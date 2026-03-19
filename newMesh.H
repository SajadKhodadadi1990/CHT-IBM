        pointField points(mesh.points());
        cellShapeList shapes(mesh.cellShapes());

        const polyBoundaryMesh& pbm = mesh.boundaryMesh();
        faceListList boundaryFaces(pbm.size()-1);
        forAll(boundaryFaces, patchi)
        {
            boundaryFaces[patchi] = pbm[patchi+1];
        }
        wordList boundaryPatchNames(pbm.names());
        wordList boundaryPatchTypes(pbm.size()-1);
        wordList boundaryPatchNames2(pbm.size()-1);
        
        forAll(boundaryPatchNames2, II)
        {
            boundaryPatchNames2[II] = boundaryPatchNames[II+1];
        }
        forAll(boundaryPatchTypes, II)
        {
            boundaryPatchTypes[II] = pbm[II+1].type();
        }
        word defaultBoundaryPatchName = "defaultFaces";
        word defaultBoundaryPatchType = emptyPolyPatch::typeName;

        fvMesh newMesh
        (
            IOobject
            (
                "newMesh",
                runTime.timeName(),
                runTime,
                Foam::IOobject::NO_READ
            ),
            Xfer<pointField>(points),
            shapes,
            boundaryFaces,
            boundaryPatchNames2,
            boundaryPatchTypes,
            defaultBoundaryPatchName,
            defaultBoundaryPatchType,
            boundaryPatchTypes,
            false
        );
        //Info<<"Mesh:"<<mesh<<endl;
       // Info<<"newMesh:"<<newMesh<<endl;
