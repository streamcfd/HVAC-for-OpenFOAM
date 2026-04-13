/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2017-2022 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "faceShading.H"
#include "fvMesh.H"
#include "boundaryRadiationProperties.H"
#include "cyclicAMIPolyPatch.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "distributedTriSurfaceMesh.H"
#include "OBJstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(faceShading, 0);
}

namespace
{
bool clipRayToBounds
(
    const Foam::treeBoundBox& globalBb,
    const Foam::point& start,
    const Foam::point& farEnd,
    const Foam::vector& clipOffset,
    Foam::point& end
)
{
    Foam::point clipPt;

    if (!globalBb.intersects(farEnd, start, clipPt))
    {
        return false;
    }

    end = clipPt + clipOffset;

    return Foam::magSqr(end - start) > Foam::magSqr(clipOffset);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::faceShading::calculate()
{
    const auto& pbm = mesh_.boundaryMesh();
    const radiation::boundaryRadiationProperties& boundaryRadiation =
        radiation::boundaryRadiationProperties::New(mesh_);

    const bitSet isOpaqueFace
    (
        selectOpaqueFaces
        (
            boundaryRadiation,
            patchIDs_,
            zoneIDs_
        )
    );

    const bitSet isPrimaryFace
    (
        selectPrimaryFaces
        (
            boundaryRadiation,
            patchIDs_,
            zoneIDs_
        )
    );

    // Find faces potentially hit by solar rays
    //  - correct normal
    //  - transmissivity < 1
    labelList hitFacesIds;
    bitSet hitFacesFlips;
    selectFaces
    (
        true,   // use normal to do first filtering
        isPrimaryFace,
        patchIDs_,
        zoneIDs_,

        hitFacesIds,
        hitFacesFlips
    );

    Info<< "Number of 'potential' direct hits : "
        << returnReduce(hitFacesIds.size(), sumOp<label>()) << endl;


    // * * * * * * * * * * * * * * *
    // Create distributedTriSurfaceMesh
    Random rndGen(653213);

    // Find potential obstructions. Include all faces that might potentially
    // block (so ignore normal)
    labelList blockingFacesIds;
    bitSet blockingFacesFlips;
    selectFaces
    (
        false,   // use normal to do first filtering
        isOpaqueFace,
        patchIDs_,
        zoneIDs_,

        blockingFacesIds,
        blockingFacesFlips
    );

    const triSurface localSurface = triangulate
    (
        blockingFacesIds,
        blockingFacesFlips
    );

    // Determine mesh bounding boxes:
    List<treeBoundBox> meshBb
    (
        1,
        treeBoundBox(mesh_.points()).extend(rndGen, 1e-3)
    );

    // Dummy bounds dictionary
    dictionary dict;
    dict.add("bounds", meshBb);
    dict.add
    (
        "distributionType",
        distributionTypeName_
    );
    dict.add("mergeDistance", SMALL);

    distributedTriSurfaceMesh surfacesMesh
    (
        IOobject
        (
            "opaqueSurface.stl",
            mesh_.time().constant(),    // directory
            "triSurface",               // instance
            mesh_.time(),               // registry
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        localSurface,
        dict
    );

    {
        autoPtr<mapDistribute> faceMap;
        autoPtr<mapDistribute> pointMap;
        surfacesMesh.distribute(meshBb, true, faceMap, pointMap);
    }

    if (debug)
    {
        surfacesMesh.searchableSurface::write();
    }

    treeBoundBox globalRayBb(boundBox(mesh_.points(), true));
    const scalar rayClipTol = max(SMALL, 1e-6*globalRayBb.mag());
    globalRayBb.grow(rayClipTol);

    const scalar maxBounding =
        returnReduce(5.0*mesh_.bounds().mag(), maxOp<scalar>());

    // Calculate index of faces which have a direct hit (local)

    // Shoot Rays
    // * * * * * * * * * * * * * * * *
    {

        DynamicField<point> start(hitFacesIds.size());
        DynamicField<point> end(start.size());
        DynamicList<label> startIndex(start.size());

        const pointField& faceCentres = mesh_.faceCentres();

        const vector d(direction_*maxBounding);
        const vector clipOffset =
            -direction_/(mag(direction_) + VSMALL)*rayClipTol;

        forAll(hitFacesIds, i)
        {
            const label facei = hitFacesIds[i];
            const point& fc = faceCentres[facei];

            const point startPt(fc + clipOffset);
            point endPt(fc - d);

            // Fall back to the original long ray if clipping fails.
            // Treating clip failures as "visible" can leak shadowed faces.
            clipRayToBounds
            (
                globalRayBb,
                startPt,
                fc - d,
                clipOffset,
                endPt
            );

            start.append(startPt);
            startIndex.append(facei);
            end.append(endPt);
        }

        List<pointIndexHit> hitInfo(startIndex.size());
        surfacesMesh.findLineAny(start, end, hitInfo);

        // Collect the rays which has 'only one not wall' obstacle between
        // start and end.
        // If the ray hit itself get stored in dRayIs

        label nVisible = 0;
        forAll(hitInfo, rayI)
        {
            if (!hitInfo[rayI].hit())
            {
                nVisible++;
            }
        }
        rayStartFaces_.setSize(nVisible);
        nVisible = 0;

        forAll(hitInfo, rayI)
        {
            if (!hitInfo[rayI].hit())
            {
                rayStartFaces_[nVisible++] = startIndex[rayI];
            }
        }

        // Plot all rays between visible faces.
        if (debug)
        {
            writeRays
            (
                mesh_.time().path()/"allVisibleFaces.obj",
                end,
                start
            );
        }

        start.clear();
        startIndex.clear();
        end.clear();
    }

    rayStartTransmissivity_.setSize(rayStartFaces_.size());
    rayStartTransmissivity_ = scalar(1);

    bitSet isTransFace(mesh_.nFaces(), false);
    scalarField faceTransmissivity(mesh_.nFaces(), scalar(1));

    for (const label patchi : patchIDs_)
    {
        const auto& pp = pbm[patchi];
        const tmp<scalarField> tt = boundaryRadiation.transmissivity(patchi);
        const scalarField& t = tt.cref();

        forAll(t, i)
        {
            const scalar ti = max(min(t[i], scalar(1)), scalar(0));
            const label meshFacei = i + pp.start();

            faceTransmissivity[meshFacei] = ti;
            isTransFace[meshFacei] = (ti > 0.0 && ti < 1.0);
        }
    }

    const auto& fzs = mesh_.faceZones();
    for (const label zonei : zoneIDs_)
    {
        const auto& fz = fzs[zonei];
        const tmp<scalarField> tt = boundaryRadiation.zoneTransmissivity
        (
            zonei,
            fz
        );
        const scalarField& t = tt.cref();

        forAll(t, i)
        {
            const scalar ti = max(min(t[i], scalar(1)), scalar(0));
            const label meshFacei = fz[i];

            faceTransmissivity[meshFacei] = ti;
            isTransFace[meshFacei] = (ti > 0.0 && ti < 1.0);
        }
    }

    labelList transmissiveFacesIds;
    bitSet transmissiveFacesFlips;
    selectFaces
    (
        false,   // use normal to do first filtering
        isTransFace,
        patchIDs_,
        zoneIDs_,

        transmissiveFacesIds,
        transmissiveFacesFlips
    );

    const label nTransmissiveFaces =
        returnReduce(transmissiveFacesIds.size(), sumOp<label>());
    const label nVisibleFaces = returnReduce(rayStartFaces_.size(), sumOp<label>());

    if (nTransmissiveFaces && nVisibleFaces)
    {
        labelList triToFace;
        const triSurface transmissiveSurface = triangulate
        (
            transmissiveFacesIds,
            transmissiveFacesFlips,
            &triToFace
        );

        distributedTriSurfaceMesh transmissiveSurfacesMesh
        (
            IOobject
            (
                "transmissiveSurface.stl",
                mesh_.time().constant(),    // directory
                "triSurface",               // instance
                mesh_.time(),               // registry
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            transmissiveSurface,
            dict
        );

        {
            autoPtr<mapDistribute> faceMap;
            autoPtr<mapDistribute> pointMap;
            transmissiveSurfacesMesh.distribute(meshBb, true, faceMap, pointMap);

            if (faceMap.valid())
            {
                faceMap().distribute(triToFace);
            }
        }

        const pointField& faceCentres = mesh_.faceCentres();
        const vector d(direction_*maxBounding);
        const vector clipOffset =
            -direction_/(mag(direction_) + VSMALL)*rayClipTol;

        DynamicField<point> start(rayStartFaces_.size());
        DynamicField<point> end(start.size());
        DynamicList<label> rayIndex(start.size());

        forAll(rayStartFaces_, i)
        {
            const label facei = rayStartFaces_[i];
            const point& fc = faceCentres[facei];

            const point startPt(fc + clipOffset);
            point endPt(fc - d);

            // Fall back to the original long ray if clipping fails so
            // transmissive layers are still tracked correctly.
            clipRayToBounds
            (
                globalRayBb,
                startPt,
                fc - d,
                clipOffset,
                endPt
            );

            start.append(startPt);
            end.append(endPt);
            rayIndex.append(i);
        }

        List<List<pointIndexHit>> allHitInfo(start.size());
        transmissiveSurfacesMesh.findLineAll(start, end, allHitInfo);

        DynamicList<pointIndexHit> flatHitInfo;
        DynamicList<label> flatRayIndex;

        forAll(allHitInfo, rayI)
        {
            for (const pointIndexHit& hit : allHitInfo[rayI])
            {
                if (hit.hit())
                {
                    flatHitInfo.append(hit);
                    flatRayIndex.append(rayI);
                }
            }
        }

        const label nFlatHits = returnReduce(flatHitInfo.size(), sumOp<label>());

        if (nFlatHits)
        {
            const List<pointIndexHit> flatHits(flatHitInfo);
            const labelList flatRays(flatRayIndex);

            labelList triangleIndex;
            autoPtr<mapDistribute> mapPtr
            (
                transmissiveSurfacesMesh.localQueries
                (
                    flatHits,
                    triangleIndex
                )
            );

            const mapDistribute& map = mapPtr();

            scalarField flatTransmissivity(triangleIndex.size(), scalar(1));
            forAll(triangleIndex, i)
            {
                const label trii = triangleIndex[i];
                if (trii >= 0)
                {
                    const label meshFacei = triToFace[trii];
                    flatTransmissivity[i] = faceTransmissivity[meshFacei];
                }
            }

            map.reverseDistribute(flatHits.size(), flatTransmissivity);

            forAll(flatTransmissivity, i)
            {
                const label rayI = rayIndex[flatRays[i]];
                rayStartTransmissivity_[rayI] *= flatTransmissivity[i];
            }
        }
    }

    if (debug)
    {
        auto thitFaces = surfaceScalarField::New
        (
            "hitFaces",
            IOobject::NO_REGISTER,
            mesh_,
            dimensionedScalar(dimless, Zero)
        );
        auto& hitFaces = thitFaces.ref();

        surfaceScalarField::Boundary& hitFacesBf = hitFaces.boundaryFieldRef();

        hitFacesBf = 0.0;
        for (const label facei : rayStartFaces_)
        {
            const label patchID = pbm.whichPatch(facei);
            if (patchID == -1)
            {
                hitFaces[facei] = 1.0;
            }
            else
            {
                const polyPatch& pp = pbm[patchID];
                hitFacesBf[patchID][facei - pp.start()] = 1.0;
            }
        }
        hitFaces.write();
    }

    Info<< "Total number of hit faces : "
        << returnReduce(rayStartFaces_.size(), sumOp<label>()) << endl;
}


Foam::triSurface Foam::faceShading::triangulate
(
    const labelUList& faceIDs,
    const bitSet& flipMap,
    labelList* triToFace
) const
{
    if (faceIDs.size() != flipMap.size())
    {
        FatalErrorInFunction << "Size problem :"
            << "faceIDs:" << faceIDs.size()
            << "flipMap:" << flipMap.size()
            << exit(FatalError);
    }

    const auto& points = mesh_.points();
    const auto& faces = mesh_.faces();
    const auto& bMesh = mesh_.boundaryMesh();
    const auto& fzs = mesh_.faceZones();

    // Patching of surface:
    // - non-processor patches
    // - faceZones
    // Note: check for faceZones on boundary? Who has priority?
    geometricSurfacePatchList surfPatches(bMesh.nNonProcessor()+fzs.size());
    labelList patchID(mesh_.nFaces(), -1);
    {
        label newPatchi = 0;
        for (label patchi = 0; patchi < bMesh.nNonProcessor(); ++patchi)
        {
            const auto& pp = bMesh[patchi];

            surfPatches[newPatchi] = geometricSurfacePatch
            (
                pp.name(),
                newPatchi,
                pp.type()
            );
            SubList<label>
            (
                patchID,
                pp.size(),
                pp.start()
            ) = newPatchi;

            newPatchi++;
        }
        for (const auto& fz : fzs)
        {
            surfPatches[newPatchi] = geometricSurfacePatch
            (
                fz.name(),
                newPatchi,
                fz.type()
            );
            UIndirectList<label>(patchID, fz) = newPatchi;

            newPatchi++;
        }
    }


    // Storage for surfaceMesh. Size estimate.
    DynamicList<labelledTri> triangles(2*faceIDs.size());
    DynamicList<label> dynTriToFace(2*faceIDs.size());
    // Work array
    faceList triFaces;

    forAll(faceIDs, i)
    {
        const label facei = faceIDs[i];
        const bool flip = flipMap[i];
        const label patchi = patchID[facei];
        const face& f = faces[facei];

        // Triangulate face
        triFaces.setSize(f.nTriangles(points));
        label nTri = 0;
        f.triangles(points, nTri, triFaces);

        for (const face& f : triFaces)
        {
            if (!flip)
            {
                triangles.append(labelledTri(f[0], f[1], f[2], patchi));
            }
            else
            {
                triangles.append(labelledTri(f[0], f[2], f[1], patchi));
            }

            if (triToFace)
            {
                dynTriToFace.append(facei);
            }
        }
    }

    triangles.shrink();

    // Create globally numbered tri surface
    triSurface rawSurface(triangles, mesh_.points());

    // Create locally numbered tri surface
    triSurface surface
    (
        rawSurface.localFaces(),
        rawSurface.localPoints()
    );

    // Add patch names to surface
    surface.patches().transfer(surfPatches);

    if (triToFace)
    {
        dynTriToFace.shrink();
        *triToFace = std::move(dynTriToFace);
    }

    return surface;
}


Foam::bitSet Foam::faceShading::selectOpaqueFaces
(
    const radiation::boundaryRadiationProperties& boundaryRadiation,
    const labelUList& patchIDs,
    const labelUList& zoneIDs
) const
{
    const auto& pbm = mesh_.boundaryMesh();

    bitSet isOpaqueFace(mesh_.nFaces(), false);

    // Check selected patches
    for (const label patchi : patchIDs)
    {
        const auto& pp = pbm[patchi];
        tmp<scalarField> tt = boundaryRadiation.transmissivity(patchi);
        const scalarField& t = tt.cref();

        forAll(t, i)
        {
            isOpaqueFace[i + pp.start()] = (t[i] == 0.0);
        }
    }

    // Check selected faceZones
    const auto& fzs = mesh_.faceZones();

    for (const label zonei : zoneIDs)
    {
        const auto& fz = fzs[zonei];

        //- Note: slice mesh face centres preferentially
        tmp<scalarField> tt = boundaryRadiation.zoneTransmissivity
        (
            zonei,
            fz
        );
        const scalarField& t = tt.cref();

        forAll(t, i)
        {
            isOpaqueFace[fz[i]] = (t[i] == 0.0);
        }
    }

    return isOpaqueFace;
}


Foam::bitSet Foam::faceShading::selectPrimaryFaces
(
    const radiation::boundaryRadiationProperties& boundaryRadiation,
    const labelUList& patchIDs,
    const labelUList& zoneIDs
) const
{
    const auto& pbm = mesh_.boundaryMesh();

    bitSet isPrimaryFace(mesh_.nFaces(), false);

    for (const label patchi : patchIDs)
    {
        const auto& pp = pbm[patchi];
        tmp<scalarField> tt = boundaryRadiation.transmissivity(patchi);
        const scalarField& t = tt.cref();

        forAll(t, i)
        {
            const scalar ti = max(min(t[i], scalar(1)), scalar(0));
            isPrimaryFace[i + pp.start()] = (ti < 1.0);
        }
    }

    const auto& fzs = mesh_.faceZones();

    for (const label zonei : zoneIDs)
    {
        const auto& fz = fzs[zonei];
        tmp<scalarField> tt = boundaryRadiation.zoneTransmissivity
        (
            zonei,
            fz
        );
        const scalarField& t = tt.cref();

        forAll(t, i)
        {
            const scalar ti = max(min(t[i], scalar(1)), scalar(0));
            isPrimaryFace[fz[i]] = (ti < 1.0);
        }
    }

    return isPrimaryFace;
}


void Foam::faceShading::selectFaces
(
    const bool useNormal,
    const bitSet& isCandidateFace,
    const labelUList& patchIDs,
    const labelUList& zoneIDs,

    labelList& faceIDs,
    bitSet& flipMap
) const
{
    const auto& pbm = mesh_.boundaryMesh();

    bitSet isSelected(mesh_.nFaces());
    DynamicList<label> dynFaces(mesh_.nBoundaryFaces());
    bitSet isFaceFlipped(mesh_.nFaces());

    // Add patches
    for (const label patchi : patchIDs)
    {
        const auto& pp = pbm[patchi];
        const vectorField& n = pp.faceNormals();

        forAll(n, i)
        {
            const label meshFacei = i + pp.start();
            if
            (
                isCandidateFace[meshFacei]
             && (
                    !useNormal
                 || ((direction_ & n[i]) > 0)
                )
            )
            {
                isSelected.set(meshFacei);
                isFaceFlipped[meshFacei] = false;
                dynFaces.append(meshFacei);
            }
        }
    }


    // Add faceZones
    const auto& fzs = mesh_.faceZones();

    for (const label zonei : zoneIDs)
    {
        const auto& fz = fzs[zonei];
        const primitiveFacePatch& pp = fz();
        const vectorField& n = pp.faceNormals();

        forAll(n, i)
        {
            const label meshFacei = fz[i];

            if
            (
                !isSelected[meshFacei]
             && isCandidateFace[meshFacei]
             && (
                    !useNormal
                 || ((direction_ & n[i]) > 0)
                )
            )
            {
                isSelected.set(meshFacei);
                dynFaces.append(meshFacei);
                isFaceFlipped[meshFacei] = fz.flipMap()[i];
            }
        }
    }
    faceIDs = std::move(dynFaces);
    flipMap = bitSet(isFaceFlipped, faceIDs);
}


void Foam::faceShading::writeRays
(
    const fileName& fName,
    const DynamicField<point>& endCf,
    const pointField& myFc
)
{
    OBJstream os(fName);

    Pout<< "Dumping rays to " << os.name() << endl;

    forAll(myFc, facei)
    {
        os.writeLine(myFc[facei], endCf[facei]);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faceShading::faceShading
(
    const fvMesh& mesh,
    const vector& dir,
    const word& distributionTypeName
)
:
    mesh_(mesh),
    patchIDs_(nonCoupledPatches(mesh)),
    zoneIDs_(0),
    direction_(dir),
    rayStartFaces_(0),
    rayStartTransmissivity_(0),
    distributionTypeName_(distributionTypeName)
{
    calculate();
}


Foam::faceShading::faceShading
(
    const fvMesh& mesh,
    const labelList& patchIDs,
    const labelList& zoneIDs,
    const vector& dir,
    const word& distributionTypeName
)
:
    mesh_(mesh),
    patchIDs_(patchIDs),
    zoneIDs_(zoneIDs),
    direction_(dir),
    rayStartFaces_(0),
    rayStartTransmissivity_(0),
    distributionTypeName_(distributionTypeName)
{
    calculate();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::faceShading::nonCoupledPatches(const polyMesh& mesh)
{
    const auto& pbm = mesh.boundaryMesh();

    DynamicList<label> ncPatches;
    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];
        if (!pp.coupled() && !isA<cyclicAMIPolyPatch>(pp))
        {
            ncPatches.append(patchi);
        }
    }
    return ncPatches;
}


void Foam::faceShading::correct()
{
    calculate();
}


// ************************************************************************* //
