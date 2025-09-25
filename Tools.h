#ifndef TOOLS_H
#define TOOLS_H

#include <deque>
#include <boost/property_tree/json_parser.hpp>
#include <vector>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include "Custom_Constrained_Delaunay_Triangulation_2.h"
#include <CGAL/draw_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/convex_hull_2.h>
#include <gmp.h>
#include <pthread.h>
#include <cmath>
using namespace std;

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef CGAL::Line_2<K> Line_2;
typedef K::Point_2                                          Point_2;
typedef K::Segment_2                                        Segment_2;
typedef vector<Point_2> Points;
typedef CGAL::Exact_predicates_tag                          Itag;
typedef Custom_Constrained_Delaunay_triangulation_2<K, CGAL::Default, Itag> CDT;
typedef CDT::Point Point;
typedef CDT::Edge  Edge;
typedef CDT::Face_handle                                    Face_handle;
typedef CDT::Vertex_handle                                  Vertex_handle;
typedef CGAL::Triangulation_vertex_base_2<K>                      Vb;
typedef CGAL::Constrained_triangulation_face_base_2<K>            Fb;
typedef CGAL::Triangle_2<K>                                 Triangle;
typedef CGAL::Angle                                     Angle;
typedef CDT::Finite_faces_iterator                      Face_handle_iterator;

struct PairInts
{
    int x, y;
    PairInts(int xi, int yi) : x(xi), y(yi) {}
};

struct Input1
{
    string instance_uid;
    int num_points;
    deque<int> points_x, points_y, region_boundary;
    int num_constraints;
    deque<PairInts> additional_constraints;
    void parse(string filename);
};

struct Parameters {
    double alpha, beta, xi, psi, lambda;
    int L, kappa;
    Parameters() {alpha = 0; beta = 0; xi = 0; psi = 0; lambda = 0; kappa = 0; L = 0;}
    Parameters(double alpha, double beta, double xi, int psi, int lambda, int kappa, int L) {
        this->alpha = alpha;
        this->beta = beta;
        this->xi = xi;
        this->psi = psi;
        this->lambda = lambda;
        this->kappa = kappa;
        this->L = L;
    }
};

struct Input2 {
    Input1 input1;
    string method, randomization;
    Parameters parameters;
    bool delauney;
    int seed;
    void parse(string filename);
};

double average_area(const CDT& cdt);
double get_rational(const K::FT& coord);
bool is_obtuse(const CDT &cdt, Face_handle fh);
bool is_finite(const CDT& cdt, Face_handle fh);
int count_obtuse_faces(const CDT& cdt);
void flip_edges(CDT &cdt);
int add_steiner_points_merge(CDT& cdt, const int, Input1& input1, vector<Point>& points);
int add_steiner_points_center(CDT& cdt, const int, Input1& input1, vector<Point>& points);
int add_steiner_points_centroid(CDT& cdt, const int, Input1& input1, vector<Point>& points);
int add_steiner_points_middle(CDT& cdt, const int, Input1& input1, vector<Point>& points);
int add_steiner_points_projection(CDT& cdt, const int, Input1& input1, vector<Point>& points);
Point get_centroid(const Polygon_2& polygon);
Point find_middle_of_largest_edge(CDT& cdt, Face_handle f);
Point find_projection_of_largest_edge(CDT& cdt, Face_handle f);
void print_rational(const K::FT& coord, ofstream& fout);
void print_result(string filename, Input2& input, CDT& cdt, bool extras = true);
int locate_index(CDT& cdt, Point &p);
bool inside_boundaries(Input1& input1, Point &p, vector<Point> points);
//------------ part 2 -------------
bool steiner_method_merge(CDT &cdt, Input1 &input1, vector<Point> &points, Face_handle_iterator f);
bool steiner_method_center(CDT &cdt, Input1 &input1, vector<Point> &points, Face_handle_iterator f);
bool steiner_method_centroid(CDT &cdt, Input1 &input1, vector<Point> &points, Face_handle_iterator f);
bool steiner_method_middle(CDT &cdt, Input1 &input1, vector<Point> &points, Face_handle_iterator f);
bool steiner_method_projection(CDT &cdt, Input1 &input1, vector<Point> &points, Face_handle_iterator f);
int LS_minimization(CDT& cdt, Input2& input, vector<Point>& points);
int simulated_annealing(CDT& cdt, Input2& input, vector<Point>& points);
double energy(CDT& cdt, Input2& input);
struct AntArgs {
    Input2& input;
    int ant_id, selector, cycle_id;
    CDT cdt_copy;
    double* ti;
    vector<Point>& points;
    AntArgs(int id, CDT& init, Input2& input2, double* t, vector<Point>& cdt_points, int cid)
        : ant_id(id), cdt_copy(init), input(input2), ti(t), points(cdt_points) { selector = -1; cycle_id = cid; }
};
int ant_colony(CDT& cdt, Input2& input, vector<Point>& points);
double find_distance_of_largest_edge(CDT& cdt, Face_handle f);
double find_height_of_largest_edge(CDT& cdt, Face_handle f);
double get_r(CDT& cdt, Face_handle f);
double heuristic_information_merge(CDT& cdt, Face_handle f);
double heuristic_information_circumcenter(CDT& cdt, Face_handle f);
double heuristic_information_centroid(CDT& cdt, Face_handle f);
double heuristic_information_middle(CDT& cdt, Face_handle f);
double heuristic_information_projection(CDT& cdt, Face_handle f);
void* ant_function(void* args);
Face_handle_iterator random_picker(CDT& cdt);
void remove_all_constraints(CDT& cdt);
vector<pair<int, int>> save_all_constraints(CDT& cdt);
void load_all_constraints(CDT& cdt, vector<pair<int, int>>& constraints, const vector<Point>& points);
vector<pair<int, int>> polygon_to_constraints(CDT& cdt, Polygon_2& polygon);
//------------ part 3 -------------
enum Category {
    A,
    B,
    C,
    D,
    E
};
Category get_CDT_Category(CDT& cdt, Input2& input, vector<Point>& points);
int add_steiner_points_random(CDT& cdt, const int, Input1& input1, vector<Point>& points);
bool steiner_method_random(CDT &cdt, Input1 &input1, vector<Point> &points, Face_handle_iterator f);
double convergence_rate();
#endif
