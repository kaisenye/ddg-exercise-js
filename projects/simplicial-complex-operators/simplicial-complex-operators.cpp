// Implement member functions for SimplicialComplexOperators class.
#include "simplicial-complex-operators.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

// Set up Macros 
typedef Eigen::SparseMatrix<size_t> SparseMat;
typedef Eigen::Triplet<size_t> Triplet;


/*
 * Assign a unique index to each vertex, edge, and face of a mesh.
 * All elements are 0-indexed.
 *
 * Input: None. Access geometry via the member variable <geometry>, and pointer to the mesh via <mesh>.
 * Returns: None.
 */
void SimplicialComplexOperators::assignElementIndices() {

    // Needed to access geometry->vertexIndices, etc. as cached quantities.
    // Not needed if you're just using v->getIndex(), etc.
    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();

    // You can set the index field of a vertex via geometry->vertexIndices[v], where v is a Vertex object (or an
    // integer). Similarly you can do edges and faces via geometry->edgeIndices, geometry->faceIndices, like so:
    size_t idx = 0;
    for (Vertex v : mesh->vertices()) {
        idx = geometry->vertexIndices[v];
    }

    for (Edge e : mesh->edges()) {
        idx = geometry->edgeIndices[e];
    }

    for (Face f : mesh->faces()) {
        idx = geometry->faceIndices[f];
    }

    // You can more easily get the indices of mesh elements using the function getIndex(), albeit less efficiently and
    // technically less safe (although you don't need to worry about it), like so:
    //
    //      v.getIndex()
    //
    // where v can be a Vertex, Edge, Face, Halfedge, etc. For example:

    for (Vertex v : mesh->vertices()) {
        idx = v.getIndex(); // == geometry->vertexIndices[v])
    }

    // Geometry Central already sets the indices for us, though, so this function is just here for demonstration.
    // You don't have to do anything :)
}

/*
 * Construct the unsigned vertex-edge adjacency matrix A0.
 *
 * Input:
 * Returns: The sparse vertex-edge adjacency matrix which gets stored in the global variable A0.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildVertexEdgeAdjacencyMatrix() const {

    // TODO
    // Note: You can build an Eigen sparse matrix from triplets, then return it as a Geometry Central SparseMatrix.
    // See <https://eigen.tuxfamily.org/dox/group__TutorialSparse.html> for documentation.
    

    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();

    std::vector<Triplet> tripletList;    //This require the macros "Triplet" declaration on line 9
    size_t numberRows = mesh->nEdges();
    size_t numberColumns = mesh->nVertices();
    SparseMatrix<size_t> matrix(numberRows,numberColumns);

    for (Edge e : mesh->edges()) {
        tripletList.push_back( 
            Triplet(geometry->edgeIndices[e], 
            geometry->vertexIndices[e.firstVertex()], 1)
        );
        tripletList.push_back(
            Triplet(geometry->edgeIndices[e], 
            geometry->vertexIndices[e.secondVertex()], 1)
        );
    }
    matrix.setFromTriplets(tripletList.begin(), tripletList.end());

    return matrix; // placeholder
}

/*
 * Construct the unsigned face-edge adjacency matrix A1.
 *
 * Input:
 * Returns: The sparse face-edge adjacency matrix which gets stored in the global variable A1.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildFaceEdgeAdjacencyMatrix() const {

    // TODO

    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();

    std::vector<Triplet> tripletList;
    size_t numberRows = mesh->nFaces();
    size_t numberColumns = mesh->nEdges();
    SparseMat matrix(numberRows, numberColumns);

    for (Face f : mesh->faces()) {
        for (Edge e : f.adjacentEdges()) {
            tripletList.push_back(
                Triplet(geometry->faceIndices[f],
                    geometry->edgeIndices[e], 1)
            );
        }        
    }
    matrix.setFromTriplets(tripletList.begin(), tripletList.end());
    return matrix; // placeholder
}

/*
 * Construct a vector encoding the vertices in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |V|, where |V| = # of vertices in the mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildVertexVector(const MeshSubset& subset) const {
    // TODO

    Vector<size_t> vect;
    vect = Vector<size_t>::Zero(mesh->nVertices());
    for (auto it = subset.vertices.begin(); it != subset.vertices.end(); it++) {  //it = iterator
        vect(*it) = 1;
    }

    return vect;
}

/*
 * Construct a vector encoding the edges in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |E|, where |E| = # of edges in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildEdgeVector(const MeshSubset& subset) const {

    // TODO

    Vector<size_t> vect;
    vect = Vector<size_t>::Zero(mesh->nEdges());
    for (auto it = subset.edges.begin(); it != subset.edges.end(); it++) {  //it = iterator
        vect(*it) = 1;
    }

    return vect;
}

/*
 * Construct a vector encoding the faces in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |F|, where |F| = # of faces in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildFaceVector(const MeshSubset& subset) const {

    // TODO
    Vector<size_t> vect;
    vect = Vector<size_t>::Zero(mesh->nFaces());
    for (auto it = subset.faces.begin(); it != subset.faces.end(); it++) {  //it = iterator
        vect(*it) = 1;
    }

    return vect;
}

/*
 * Compute the simplicial star St(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The star of the given subset.
 */
MeshSubset SimplicialComplexOperators::star(const MeshSubset& subset) const {

    // TODO

    Vector<size_t> v_vect, e_vect, extra_edge, extra_face;
    MeshSubset star;
    star = subset;

    v_vect = buildVertexVector(star);
    extra_edge=A0*v_vect;
    for (size_t i = 0; i < mesh->nEdges(); i++) {
        if (extra_edge(i) != 0) {
            star.addEdge(i);
        }
    }

    e_vect = buildEdgeVector(star);
    extra_face = A1*e_vect;
    for (size_t i = 0; i < mesh->nFaces(); i++) {
        if (extra_face(i) != 0) {
            star.addFace(i);
        }
    }

    return star; // placeholder
}


/*
 * Compute the closure Cl(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The closure of the given subset.
 */
MeshSubset SimplicialComplexOperators::closure(const MeshSubset& subset) const {

    // TODO
    MeshSubset closure;
    closure = subset;



    return closure; // placeholder
}

/*
 * Compute the link Lk(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The link of the given subset.
 */
MeshSubset SimplicialComplexOperators::link(const MeshSubset& subset) const {

    // TODO
    return subset; // placeholder
}

/*
 * Return true if the selected subset is a simplicial complex, false otherwise.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: True if given subset is a simplicial complex, false otherwise.
 */
bool SimplicialComplexOperators::isComplex(const MeshSubset& subset) const {

    // TODO
    return false; // placeholder
}

/*
 * Check if the given subset S is a pure simplicial complex. If so, return the degree of the complex. Otherwise, return
 * -1.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: int representing the degree of the given complex (-1 if not pure)
 */
int SimplicialComplexOperators::isPureComplex(const MeshSubset& subset) const {

    // TODO
    return -1; // placeholder
}

/*
 * Compute the set of simplices contained in the boundary bd(S) of the selected subset S of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The boundary of the given subset.
 */
MeshSubset SimplicialComplexOperators::boundary(const MeshSubset& subset) const {

    // TODO
    return subset; // placeholder
}