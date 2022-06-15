"use strict";

/**
 * This class contains methods to build common {@link https://cs.cmu.edu/~kmcrane/Projects/DDG/paper.pdf discrete exterior calculus} operators.
 * @memberof module:Core
 */
class DEC {
	/**
	 * Builds a sparse diagonal matrix encoding the Hodge operator on 0-forms.
	 * By convention, the area of a vertex is 1.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} vertexIndex A dictionary mapping each vertex of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildHodgeStar0Form(geometry, vertexIndex) {
		// TODO
		
		return SparseMatrix.identity(1, 1); // placeholder
	}

	/**
	 * Builds a sparse diagonal matrix encoding the Hodge operator on 1-forms.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} edgeIndex A dictionary mapping each edge of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildHodgeStar1Form(geometry, edgeIndex) {
		// TODO

		return SparseMatrix.identity(1, 1); // placeholder
	}

	/**
	 * Builds a sparse diagonal matrix encoding the Hodge operator on 2-forms.
	 * By convention, the area of a vertex is 1.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} faceIndex A dictionary mapping each face of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildHodgeStar2Form(geometry, faceIndex) {
		// TODO

		return SparseMatrix.identity(1, 1); // placeholder
	}

	/**
	 * Builds a sparse matrix encoding the exterior derivative on 0-forms.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} edgeIndex A dictionary mapping each edge of a mesh to a unique index.
	 * @param {Object} vertexIndex A dictionary mapping each vertex of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildExteriorDerivative0Form(geometry, edgeIndex, vertexIndex) {
		// TODO

		// get matrix rows and columns number
		let numnberRows = Object.keys(edgeIndex).length;
		let numberColumns = Object.keys(vertexIndex).length;

		// build a tripletList filled with 0s
		let tripletList = new Triplet(numnberRows, numberColumns);

		// main loop
		for (let e of geometry.mesh.edges) {
			tripletList.addEntry(
					1, edgeIndex[e], 
					vertexIndex[e.halfedge.vertex]
			);
			tripletList.addEntry(
					-1, edgeIndex[e],
					vertexIndex[e.halfedge.twin.vertex]
			);
		}
		let matrix = SparseMatrix.fromTriplet(tripletList);

		return matrix; // placeholder
	}

	/**
	 * Builds a sparse matrix encoding the exterior derivative on 1-forms.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} faceIndex A dictionary mapping each face of a mesh to a unique index.
	 * @param {Object} edgeIndex A dictionary mapping each edge of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildExteriorDerivative1Form(geometry, faceIndex, edgeIndex) {
		// TODO
		// get matrix rows and columns number
		let numnberRows = Object.keys(faceIndex).length;
		let numberColumns = Object.keys(edgeIndex).length;

		// build a triplet list
		let tripletList = new Triplet(numnberRows, numberColumns);

		// main loop
		for (let f of geometry.mesh.faces) {
			for (let he of f.adjacentHalfedges()) {
				if (he.vertex != he.edge.halfedge.vertex) {
					tripletList.addEntry(
						-1, faceIndex[f],
						edgeIndex[he.edge]
					)
				} else {
					tripletList.addEntry(
						1, faceIndex[f],
						edgeIndex[he.edge]
					)
				}
			}
		}


		let matrix = SparseMatrix.fromTriplet(tripletList);

		return matrix;
	}
}
