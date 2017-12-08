//=============================================================================
//
//   Code framework for the lecture
//
//   "Digital 3D Geometry Processing"
//
//   Copyright (C) 2017 by Computer Graphics and Geometry Laboratory,
//         EPF Lausanne
//
//-----------------------------------------------------------------------------
#define _USE_MATH_DEFINES
#include "mesh_processing.h"
#include <cmath>
#include <set>
#include <map>

namespace mesh_processing {

using surface_mesh::Point;
using surface_mesh::Scalar;
using surface_mesh::Color;
using std::min;
using std::max;
using std::cout;
using std::endl;

MeshProcessing::MeshProcessing(const string& filename) {
    load_mesh(filename);
}

void MeshProcessing::deformation() {
	deformation_axis(0);
	deformation_axis(1);
	deformation_axis(2);
}

void MeshProcessing::deformation_axis(int mode) {
	const int N = mesh_.n_vertices();

	calc_weights();
	auto cotan = mesh_.edge_property<Scalar>("e:weight");

	Eigen::SparseMatrix<double> L(N, N);
	Eigen::MatrixXd rhs(Eigen::MatrixXd::Zero(N, 1));
	std::vector< Eigen::Triplet<double> > triplets_L;

	
	// ! Adding itearators
	Mesh::Halfedge_around_vertex_circulator vh_c, vh_end;
	Mesh::Vertex neighbor_v;
	Mesh::Edge e;

	for (int i = 0; i < N; ++i) {

		// ------------- IMPLEMENT HERE ---------
		// Set up Laplace-Beltrami matrix of the mesh
		// ------------- IMPLEMENT HERE ---------	

		Mesh::Vertex v(i);
		auto currentEdgeAngle = 0;
		auto sumOfEdgeWeights = 0;


		//Loop through adjacent halfedges to get adjacent vertices
		vh_c = mesh_.halfedges(v);
		vh_end = vh_c;
		do {
			neighbor_v = mesh_.to_vertex(*vh_c);
			e = mesh_.find_edge(v, neighbor_v);
			currentEdgeAngle = cotan[e];
			sumOfEdgeWeights += currentEdgeAngle;
		
		} while (++vh_c != vh_end);

		triplets_L.push_back(Eigen::Triplet<double>(v.idx(), neighbor_v.idx(), currentEdgeAngle));
		triplets_L.push_back(Eigen::Triplet<double>(v.idx(), v.idx(), sumOfEdgeWeights));
	

		// ! Adding restrictions to B matrix
		for (std::vector<size_t>::iterator it = shifted_faces_points_indices_.begin();it != shifted_faces_points_indices_.end();++it)
		{
			if (v.idx() == *it)
			{	
				rhs(v.idx()) = displacement_[mode];
				}

			}
		
	}

	// ! Setting the L Matrix
	L.setFromTriplets(triplets_L.begin(), triplets_L.end());
	auto L_squared=L*L;
	// ------------- IMPLEMENT HERE ---------
	// Compute squared Laplace-Beltrami matrix 
	// For each fixed and shifted vertex replace the corresponding row of squared Laplace-Beltrami matrix with the constraint row
	//    Hint: to iterate through sparse matrix use Eigen::SparseMatrix<double>::InnerIterator.
	//          note that sparse matrix is stored column-wise
	//          since squared Laplace-Beltrami matrix is symmetric, 
	//          you can traverse the matrix column-wise and transpose the result
	// Solve the linear system L2x = b
	// Displace the vertices by x
	//    Hint: the input parameter of the function "mode" indicates the current displacement axis.
	//    the function deformation_axis() is called 3 times, with mode = 0, mode = 1 and mode = 3 for the axis X, Y and Z
	// ------------- IMPLEMENT HERE ---------
	fixed_faces_points_indices_;
	




	//Solving the linear system
	Eigen::SparseLU< Eigen::SparseMatrix<double> > solver(L_squared);
	if (solver.info() != Eigen::Success)
		printf("linear solver init failed.\n");

	Eigen::MatrixXd X = solver.solve(rhs);
	if (solver.info() != Eigen::Success)
		printf("linear solver failed.\n");


	//Setting the values
	Mesh::Vertex_property<Scalar> v_harmonic_function = mesh_.vertex_property<Scalar>("Harmonic", 0.0f);

	for (int i = 0; i < N; ++i)
		v_harmonic_function[Mesh::Vertex(i)] = X(i, mode);


	// clean-up
	mesh_.remove_edge_property(cotan);
}

void MeshProcessing::calc_uniform_mean_curvature() {
    Mesh::Vertex_property<Scalar> v_unicurvature =
            mesh_.vertex_property<Scalar>("v:unicurvature", 0.0f);
    Mesh::Vertex_around_vertex_circulator   vv_c, vv_end;
    Point             laplace(0.0);

    for (auto v: mesh_.vertices()) {
        Scalar curv = 0;

        if (!mesh_.is_boundary(v)) {
            laplace = Point(0.0f);
            double n = 0;
            vv_c = mesh_.vertices(v);
            vv_end = vv_c;

            do {
                laplace += (mesh_.position(*vv_c) - mesh_.position(v));
                ++n;
            } while(++vv_c != vv_end);

            laplace /= n;

            curv = 0.5f * norm(laplace);
        }
        v_unicurvature[v] = curv;
    }
}

void MeshProcessing::calc_mean_curvature() {
    Mesh::Vertex_property<Scalar>  v_curvature =
            mesh_.vertex_property<Scalar>("v:curvature", 0.0f);
    Mesh::Edge_property<Scalar> e_weight =
            mesh_.edge_property<Scalar>("e:weight", 0.0f);
    Mesh::Vertex_property<Scalar>  v_weight =
            mesh_.vertex_property<Scalar>("v:weight", 0.0f);

    Mesh::Halfedge_around_vertex_circulator vh_c, vh_end;
    Mesh::Vertex neighbor_v;
    Mesh::Edge e;
    Point laplace(0.0f, 0.0f, 0.0f);

    for (auto v: mesh_.vertices()) {
        Scalar curv = 0.0f;

        if (!mesh_.is_boundary(v)) {
            laplace = Point(0.0f, 0.0f, 0.0f);

            vh_c = mesh_.halfedges(v);
            vh_end = vh_c;

            do {
                e = mesh_.edge(*vh_c);
                neighbor_v = mesh_.to_vertex(*vh_c);
                laplace += e_weight[e] * (mesh_.position(neighbor_v) -
                                          mesh_.position(v));

            } while(++vh_c != vh_end);

            laplace *= v_weight[v];
            curv = 0.5f * norm(laplace);
        }
        v_curvature[v] = curv;
    }
}

void MeshProcessing::calc_gauss_curvature() {
    Mesh::Vertex_property<Scalar> v_gauss_curvature =
            mesh_.vertex_property<Scalar>("v:gauss_curvature", 0.0f);
    Mesh::Vertex_property<Scalar> v_weight =
            mesh_.vertex_property<Scalar>("v:weight", 0.0f);
    Mesh::Vertex_around_vertex_circulator vv_c, vv_c2, vv_end;
    Point d0, d1;
    Scalar angles, cos_angle;
    Scalar lb(-1.0f), ub(1.0f);

    // compute for all non-boundary vertices
    for (auto v: mesh_.vertices()) {
        Scalar curv = 0.0f;

        if (!mesh_.is_boundary(v)) {
            angles = 0.0f;

            vv_c = mesh_.vertices(v);
            vv_end = vv_c;

            do {
                vv_c2 = vv_c;
                ++ vv_c2;
                d0 = normalize(mesh_.position(*vv_c) - mesh_.position(v));
                d1 = normalize(mesh_.position(*vv_c2) - mesh_.position(v));
                cos_angle = max(lb, min(ub, dot(d0, d1)));
                angles += acos(cos_angle);
            } while(++vv_c != vv_end);

            curv = (2 * (Scalar)M_PI - angles) * 2.0f * v_weight[v];
        }
        v_gauss_curvature[v] = curv;
    }
}

void MeshProcessing::calc_weights() {
    calc_edges_weights();
    calc_vertices_weights();
}

void MeshProcessing::calc_edges_weights() {
    auto e_weight = mesh_.edge_property<Scalar>("e:weight", 0.0f);
    auto points = mesh_.vertex_property<Point>("v:point");

    Mesh::Halfedge h0, h1, h2;
    Point p0, p1, p2, d0, d1;

    for (auto e: mesh_.edges())
    {
        e_weight[e] = 0.0;

        h0 = mesh_.halfedge(e, 0);
        p0 = points[mesh_.to_vertex(h0)];

        h1 = mesh_.halfedge(e, 1);
        p1 = points[mesh_.to_vertex(h1)];

        if (!mesh_.is_boundary(h0))
        {
            h2 = mesh_.next_halfedge(h0);
            p2 = points[mesh_.to_vertex(h2)];
            d0 = p0 - p2;
            d1 = p1 - p2;
            e_weight[e] += dot(d0,d1) / norm(cross(d0,d1));
        }

        if (!mesh_.is_boundary(h1))
        {
            h2 = mesh_.next_halfedge(h1);
            p2 = points[mesh_.to_vertex(h2)];
            d0 = p0 - p2;
            d1 = p1 - p2;
            e_weight[e] += dot(d0,d1) / norm(cross(d0,d1));
        }
    }
}

void MeshProcessing::calc_vertices_weights() {
    Mesh::Face_around_vertex_circulator vf_c, vf_end;
    Mesh::Vertex_around_face_circulator fv_c;
    Scalar area;
    auto v_weight = mesh_.vertex_property<Scalar>("v:weight", 0.0f);

    for (auto v: mesh_.vertices()) {
        area = 0.0;
        vf_c = mesh_.faces(v);

        if(!vf_c) {
            continue;
        }

        vf_end = vf_c;

        do {
            fv_c = mesh_.vertices(*vf_c);

            const Point& P = mesh_.position(*fv_c);  ++fv_c;
            const Point& Q = mesh_.position(*fv_c);  ++fv_c;
            const Point& R = mesh_.position(*fv_c);

            area += norm(cross(Q-P, R-P)) * 0.5f * 0.3333f;

        } while(++vf_c != vf_end);

        v_weight[v] = 0.5 / area;
    }
}

void MeshProcessing::load_mesh(const string &filename) {
    if (!mesh_.read(filename)) {
        std::cerr << "Mesh not found, exiting." << std::endl;
        exit(-1);
    }

    cout << "Mesh "<< filename << " loaded." << endl;
    cout << "# of vertices : " << mesh_.n_vertices() << endl;
    cout << "# of faces : " << mesh_.n_faces() << endl;
    cout << "# of edges : " << mesh_.n_edges() << endl;

    // Compute the center of the mesh
    mesh_center_ = Point(0.0f, 0.0f, 0.0f);
    for (auto v: mesh_.vertices()) {
        mesh_center_ += mesh_.position(v);
    }
    mesh_center_ /= mesh_.n_vertices();

    // Compute the maximum distance from all points in the mesh and the center
    dist_max_ = 0.0f;
    for (auto v: mesh_.vertices()) {
        if (distance(mesh_center_, mesh_.position(v)) > dist_max_) {
            dist_max_ = distance(mesh_center_, mesh_.position(v));
        }
    }

    compute_mesh_properties();

    // Store the original mesh, this might be useful for some computations
    mesh_init_ = mesh_;
}

void MeshProcessing::compute_mesh_properties() {
    Mesh::Vertex_property<Point> vertex_normal =
            mesh_.vertex_property<Point>("v:normal");
    mesh_.update_face_normals();
    mesh_.update_vertex_normals();
    Mesh::Vertex_property<Color> v_color_valence =
            mesh_.vertex_property<Color>("v:color_valence",
                                         Color(1.0f, 1.0f, 1.0f));
    Mesh::Vertex_property<Color> v_color_unicurvature =
            mesh_.vertex_property<Color>("v:color_unicurvature",
                                         Color(1.0f, 1.0f, 1.0f));
    Mesh::Vertex_property<Color> v_color_curvature =
            mesh_.vertex_property<Color>("v:color_curvature",
                                         Color(1.0f, 1.0f, 1.0f));
    Mesh::Vertex_property<Color> v_color_gaussian_curv =
            mesh_.vertex_property<Color>("v:color_gaussian_curv",
                                         Color(1.0f, 1.0f, 1.0f));
	Mesh::Vertex_property<Color> v_color_laplacian =
		mesh_.vertex_property<Color>("v:color_laplacian",
			Color(1.0f, 1.0f, 1.0f));

    Mesh::Vertex_property<Scalar> vertex_valence =
            mesh_.vertex_property<Scalar>("v:valence", 0.0f);
    for (auto v: mesh_.vertices()) {
        vertex_valence[v] = mesh_.valence(v);
    }

    Mesh::Vertex_property<Scalar> v_unicurvature =
            mesh_.vertex_property<Scalar>("v:unicurvature", 0.0f);
    Mesh::Vertex_property<Scalar> v_curvature =
            mesh_.vertex_property<Scalar>("v:curvature", 0.0f);
    Mesh::Vertex_property<Scalar> v_gauss_curvature =
            mesh_.vertex_property<Scalar>("v:gauss_curvature", 0.0f);
	Mesh::Vertex_property<Scalar> v_harmonic_function =
		mesh_.vertex_property<Scalar>("v:harmonic_function_0", 0.0f);

    calc_weights();
    calc_uniform_mean_curvature();
    calc_mean_curvature();
    calc_gauss_curvature();
    color_coding(vertex_valence, &mesh_, v_color_valence, 100 /* bound */);
    color_coding(v_unicurvature, &mesh_, v_color_unicurvature);
    color_coding(v_curvature, &mesh_, v_color_curvature);
    color_coding(v_gauss_curvature, &mesh_, v_color_gaussian_curv);
	color_coding(v_harmonic_function, &mesh_, v_color_laplacian);

    // get the mesh attributes and upload them to the GPU
    int j = 0;
    unsigned int n_vertices(mesh_.n_vertices());

    // Create big matrices to send the data to the GPU with the required
    // format
    color_valence_ = Eigen::MatrixXf(3, n_vertices);
    color_unicurvature_ = Eigen::MatrixXf(3, n_vertices);
    color_curvature_ = Eigen::MatrixXf(3, n_vertices);
    color_gaussian_curv_ = Eigen::MatrixXf(3, n_vertices);
	color_laplacian_ = Eigen::MatrixXf(3, n_vertices);
    normals_ = Eigen::MatrixXf(3, n_vertices);
    points_ = Eigen::MatrixXf(3, n_vertices);	
	selection_ = Eigen::MatrixXf(3, 4);
    indices_ = MatrixXu(3, mesh_.n_faces());

    for(auto f: mesh_.faces()) {
        std::vector<float> vv(3);
        int k = 0;
        for (auto v: mesh_.vertices(f)) {
            vv[k] = v.idx();
            ++k;
        }
        indices_.col(j) << vv[0], vv[1], vv[2];
        ++j;
    }

    j = 0;
    for (auto v: mesh_.vertices()) {
        points_.col(j) << mesh_.position(v).x,
                          mesh_.position(v).y,
                          mesh_.position(v).z;

        normals_.col(j) << vertex_normal[v].x,
                           vertex_normal[v].y,
                           vertex_normal[v].z;

        color_valence_.col(j) << v_color_valence[v].x,
                                 v_color_valence[v].y,
                                 v_color_valence[v].z;

        color_unicurvature_.col(j) << v_color_unicurvature[v].x,
                                      v_color_unicurvature[v].y,
                                      v_color_unicurvature[v].z;

        color_curvature_.col(j) << v_color_curvature[v].x,
                                   v_color_curvature[v].y,
                                   v_color_curvature[v].z;

        color_gaussian_curv_.col(j) << v_color_gaussian_curv[v].x,
                                       v_color_gaussian_curv[v].y,
                                       v_color_gaussian_curv[v].z;

		color_laplacian_.col(j) << v_color_laplacian[v].x,
								   v_color_laplacian[v].y,
								   v_color_laplacian[v].z;
        ++j;
    }
}

void MeshProcessing::color_coding(Mesh::Vertex_property<Scalar> prop, Mesh *mesh,
                  Mesh::Vertex_property<Color> color_prop, int bound) {
    // Get the value array
    std::vector<Scalar> values = prop.vector();

    // discard upper and lower bound
    unsigned int n = values.size()-1;
    unsigned int i = n / bound;
    std::sort(values.begin(), values.end());
    Scalar min_value = values[i], max_value = values[n-1-i];

    // map values to colors
    for (auto v: mesh->vertices())
    {
        set_color(v, value_to_color(prop[v], min_value, max_value), color_prop);
    }
}

void MeshProcessing::set_color(Mesh::Vertex v, const Color& col,
               Mesh::Vertex_property<Color> color_prop)
{
    color_prop[v] = col;
}

Color MeshProcessing::value_to_color(Scalar value, Scalar min_value, Scalar max_value) {
    Scalar v0, v1, v2, v3, v4;
    v0 = min_value + 0.0/4.0 * (max_value - min_value);
    v1 = min_value + 1.0/4.0 * (max_value - min_value);
    v2 = min_value + 2.0/4.0 * (max_value - min_value);
    v3 = min_value + 3.0/4.0 * (max_value - min_value);
    v4 = min_value + 4.0/4.0 * (max_value - min_value);

    Color col(1.0f, 1.0f, 1.0f);

    if (value < v0) {
        col = Color(0, 0, 1);
    } else if (value > v4) {
        col = Color(1, 0, 0);
    } else if (value <= v2) {
        if (value <= v1) { // [v0, v1]
            Scalar u =  (value - v0) / (v1 - v0);
            col = Color(0, u, 1);
        } else { // ]v1, v2]
            Scalar u = (value - v1) / (v2 - v1);
            col = Color(0, 1, 1-u);
        }
    } else {
        if (value <= v3) { // ]v2, v3]
            Scalar u = (value - v2) / (v3 - v2);
            col = Color(u, 1, 0);
        } else { // ]v3, v4]
            Scalar u = (value - v3) / (v4 - v3);
            col = Color(1, 1-u, 0);
        }
    }
    return col;
}

Eigen::Vector3f MeshProcessing::get_closest_vertex(const Eigen::Vector3f & origin, const Eigen::Vector3f & direction, size_t & closest_index) {
	float min_distance = std::numeric_limits<float>::max();
	Eigen::Vector3f closest_vertex;

	for (int i = 0; i <  mesh_.n_vertices(); ++i) {
		Mesh::Vertex v(i);
		Eigen::Vector3f point;
		point << mesh_.position(v).x, mesh_.position(v).y, mesh_.position(v).z;
		float projection_length = (point - origin).dot(direction);
		Eigen::Vector3f difference = point - (origin + projection_length * direction);
		float distance = difference.norm();
		if (distance < min_distance) {
			closest_index = i;
			min_distance = distance;
			closest_vertex = point;
		}
	}
	return closest_vertex;
}

bool ray_triangle_intersection(const Eigen::Vector3f & p0, const Eigen::Vector3f & p1, const Eigen::Vector3f & p2, const Eigen::Vector3f & o, const Eigen::Vector3f & d, Eigen::Vector3f & i) {
	float epsilon = 0.0000001;
	Eigen::Vector3f e1 = p1 - p0;
	Eigen::Vector3f e2 = p2 - p0;
	Eigen::Vector3f q = d.cross(e2);
	float a = e1.dot(q);

	// the vector is parallel to the plane(the intersection is at infinity)
	if (a > -epsilon && a < epsilon) {
		return false;
	}
	float f = 1 / a;
	Eigen::Vector3f s = o - p0;
	float u = f * s.dot(q);

	// the intersection is outside of the triangle
	if (u < 0.0) {
		return false;
	}
	Eigen::Vector3f r = s.cross(e1);
	float v = f * d.dot(r);

	// the intersection is outside of the triangle
	if (v <0.0 || u + v > 1.0) {
		return false;
	}
	float t = f * e2.dot(r);
	i = o + t * d;
	return true;
}

int MeshProcessing::get_closest_face(const Eigen::Vector3f & origin, const Eigen::Vector3f & direction, int mode) {
	float min_distance = std::numeric_limits<float>::max();
	int closest_index = -1;

	//for (Mesh::Face_iterator f_it = mesh_.faces_begin(); f_it != mesh_.faces_end(); ++f_it) {
	for (size_t i = 0; i < mesh_.n_faces(); i++) {
		Mesh::Vertex_around_face_circulator fv_c = mesh_.vertices(Mesh::Face(i));
		Mesh::Vertex_around_face_circulator  fv_end = fv_c;
		std::vector<Point> points;
		do {
			Mesh::Vertex v = (*fv_c);
			Point p = mesh_.position(v);
			points.push_back(p);
		} while (++fv_c != fv_end);
		Eigen::Vector3f p0 = Eigen::Vector3f(points[0][0], points[0][1], points[0][2]);
		Eigen::Vector3f p1 = Eigen::Vector3f(points[1][0], points[1][1], points[1][2]);
		Eigen::Vector3f p2 = Eigen::Vector3f(points[2][0], points[2][1], points[2][2]);
		Eigen::Vector3f intersection;
		bool is_intersecting = ray_triangle_intersection(p0, p1, p2, origin, direction, intersection);
		if (!is_intersecting) continue;
		float distance = (origin - intersection).norm();
		if (distance < min_distance) {
			closest_index = i;
			min_distance = distance;
		}
	}
	if (closest_index >= 0) {
		Mesh::Vertex_around_face_circulator fv_c = mesh_.vertices(Mesh::Face(closest_index));
		Mesh::Vertex_around_face_circulator  fv_end = fv_c;		
		do {
			Mesh::Vertex v = (*fv_c);
			Point p = mesh_.position(v);
			if (mode == 0) {
				fixed_faces_points_.push_back(p);
				fixed_faces_points_indices_.push_back(v.idx());
			}
			if (mode == 1) {
				shifted_faces_points_.push_back(p);
				shifted_faces_points_indices_.push_back(v.idx());
			}
		} while (++fv_c != fv_end);		
	}
	return closest_index;
}

MeshProcessing::~MeshProcessing() {}
}
