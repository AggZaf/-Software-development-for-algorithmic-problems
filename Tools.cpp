#include "Tools.h"

vector<pair<int, int>> convergence; // pairs of <steiner_points, obtuse_faces>

void Input1::parse(string filename)
{
    cout << "Reading json file: " << filename << endl;
    // Short alias for this namespace
    namespace pt = boost::property_tree;

    // Create a root
    pt::ptree root;

    // Load the json file in this ptree
    pt::read_json(filename, root);

    // Read values
    num_points = root.get<int>("num_points", 0);
    instance_uid = root.get<string>("instance_uid", "");
    for (pt::ptree::value_type &point : root.get_child("points_x"))
        points_x.push_back(atoi(point.second.data().c_str()));
    for (pt::ptree::value_type &point : root.get_child("points_y"))
        points_y.push_back(atoi(point.second.data().c_str()));
    for (pt::ptree::value_type &point : root.get_child("region_boundary"))
        region_boundary.push_back(atoi(point.second.data().c_str()));
    num_constraints = root.get<int>("num_constraints", 0);
    int add_constraints[num_constraints][2];
    int x = 0;
    for (pt::ptree::value_type &row : root.get_child("additional_constraints"))
    {
        int y = 0;
        for (pt::ptree::value_type &cell : row.second)
        {
            add_constraints[x][y] = cell.second.get_value<int>();
            y++;
        }
        PairInts pi(add_constraints[x][0], add_constraints[x][1]);
        additional_constraints.push_back(pi);
        x++;
    }
    cout << "Completed!" << endl;
}

void Input2::parse(string filename)
{
    cout << "Reading json file (for part 2): " << filename << endl;
    // Short alias for this namespace
    namespace pt = boost::property_tree;

    // Create a root
    pt::ptree root;

    // Load the json file in this ptree
    pt::read_json(filename, root);

    // Read values
    input1.num_points = root.get<int>("num_points", 0);
    input1.instance_uid = root.get<string>("instance_uid", "");
    for (pt::ptree::value_type &point : root.get_child("points_x"))
        input1.points_x.push_back(atoi(point.second.data().c_str()));
    for (pt::ptree::value_type &point : root.get_child("points_y"))
        input1.points_y.push_back(atoi(point.second.data().c_str()));
    for (pt::ptree::value_type &point : root.get_child("region_boundary"))
        input1.region_boundary.push_back(atoi(point.second.data().c_str()));
    input1.num_constraints = root.get<int>("num_constraints", 0);
    int add_constraints[input1.num_constraints][2];
    int x = 0;
    for (pt::ptree::value_type &row : root.get_child("additional_constraints"))
    {
        int y = 0;
        for (pt::ptree::value_type &cell : row.second)
        {
            add_constraints[x][y] = cell.second.get_value<int>();
            y++;
        }
        PairInts pi(add_constraints[x][0], add_constraints[x][1]);
        input1.additional_constraints.push_back(pi);
        x++;
    }
    method = root.get<string>("method", "");
    string delaunay = root.get<string>("delaunay", "");
    this->delauney = false;
    if (delaunay == "true")
        this->delauney = true;

    // Iterator over all parameters
    if (root.count("parameters") == 0)
        return;
    for (pt::ptree::value_type &parameter : root.get_child("parameters"))
    {
        // Parameter is a pair of a string and a double

        // Get the label of the node
        string name = parameter.first;
        // Get the content of the node
        string value = parameter.second.data();
        if (name == "alpha")
            parameters.alpha = atof(value.c_str());
        else if (name == "beta")
            parameters.beta = atof(value.c_str());
        else if (name == "xi")
            parameters.xi = atof(value.c_str());
        else if (name == "psi")
            parameters.psi = atof(value.c_str());
        else if (name == "lambda")
            parameters.lambda = atof(value.c_str());
        else if (name == "kappa")
            parameters.kappa = atoi(value.c_str());
        else if (name == "L")
            parameters.L = atoi(value.c_str());
    }
    randomization = "false";
    cout << "Completed!" << endl;
}

bool steiner_method_merge(CDT &cdt, Input1 &input1, vector<Point> &points, Face_handle_iterator f) {
    vector<int> facesToMerge;

    //should have at least 1 obtuse neighbors
    if (is_finite(cdt, f->neighbor(0)) && !is_obtuse(cdt, f->neighbor(0)) &&
        is_finite(cdt, f->neighbor(1)) && !is_obtuse(cdt, f->neighbor(1)) &&
        is_finite(cdt, f->neighbor(2)) && !is_obtuse(cdt, f->neighbor(2)))
        return true;

    for (int n = 0; n < 3; n++) {
        Face_handle ni = f->neighbor(n);
        if (is_finite(cdt, ni) == false)
            continue;
        if (ni->is_constrained(n))
            continue;
        facesToMerge.push_back(n);
    }
    if (facesToMerge.size() == 0)
        return true;

    Polygon_2 polygon;
    // cout << "Polygon:" << endl;
    for (int i = 0; i < 3; i++) {
        const Point& p = f->vertex(i)->point();
        polygon.push_back(p);
        // cout << p << endl;
    }
    for (int face : facesToMerge) {
        Face_handle ni = f->neighbor(face);
        for (int i = 0; i < 3; i++) {
            if (f->has_vertex(ni->vertex(i)))
                continue;
            polygon.push_back(ni->vertex(i)->point());
            // cout << ni->vertex(i)->point() << endl;
        }
    }
    if (polygon.vertices().size() <= 3)
        return true;

    Point centroid = get_centroid(polygon);
    bool point_exists = false;
    for (CDT::Vertex_handle v : cdt.all_vertex_handles())
        if (v->point() == centroid) {
            point_exists = true;
            break;
        }
    if (point_exists)
        return true;
    if (!inside_boundaries(input1, centroid, points))
        return true;

    //vector<pair<int, int>> constraints_current = save_all_constraints(cdt);
    //vector<pair<int, int>> constraints_polygon = polygon_to_constraints(cdt, polygon);
    bool no_flip = false;
    CDT cdt_copy(cdt);
    // draw(cdt_copy);
    //load_all_constraints(cdt, constraints_polygon, points);
    // draw(cdt_copy);
    cdt_copy.insert(centroid);
    // draw(cdt_copy);
    //remove_all_constraints(cdt_copy);
    //load_all_constraints(cdt_copy, constraints_current, points);
    // draw(cdt_copy);

    if (count_obtuse_faces(cdt_copy) > count_obtuse_faces(cdt)) {
        cdt_copy = cdt;
        //load_all_constraints(cdt, constraints_polygon, points);
        cdt_copy.insert_no_flip(centroid);
        //remove_all_constraints(cdt_copy);
        //load_all_constraints(cdt, constraints_current, points);
        if (count_obtuse_faces(cdt_copy) > count_obtuse_faces(cdt))
            return true;
        no_flip = true;
    }
    if (no_flip) {
        //load_all_constraints(cdt, constraints_polygon, points);
        cdt.insert_no_flip(centroid);
        //remove_all_constraints(cdt);
        //load_all_constraints(cdt, constraints_current, points);
    }
    else {
        //load_all_constraints(cdt, constraints_polygon, points);
        cdt.insert(centroid);
        //remove_all_constraints(cdt);
        //load_all_constraints(cdt, constraints_current, points);
    }
    return false;
}

int add_steiner_points_merge(CDT& cdt, const int max_limit, Input1& input1, vector<Point>& points) {
    int steiner_points = 0;
    for (Face_handle_iterator f = cdt.finite_faces_begin(); f != cdt.finite_faces_end(); f++)
    {
        if (is_obtuse(cdt, f)) {
            if (steiner_method_merge(cdt, input1, points, f))
                continue;
            steiner_points++;
            convergence.push_back(make_pair(cdt.finite_vertex_handles().size() - points.size(), count_obtuse_faces(cdt)));
            //cout << "Inserted (merge): " << centroid << endl;
            //cout << "Average area per triangle: " << average_area(cdt) << endl;
            if (steiner_points > max_limit)
                return steiner_points;
            f = cdt.finite_faces_begin();
            for (int i = 0; i < rand() % cdt.number_of_faces() - 4; i++)
                f++;
        }
    }
    return steiner_points;
}

bool steiner_method_center(CDT &cdt, Input1 &input1, vector<Point> &points, Face_handle_iterator f) {
    Triangle t = cdt.triangle(f);
    Point p1 = f->vertex(0)->point();
    Point p2 = f->vertex(1)->point();
    Point p3 = f->vertex(2)->point();

    Point new_point = CGAL::circumcenter(t);
    if (!inside_boundaries(input1, new_point, points)) {
        new_point = CGAL::centroid(t);
    }

    bool point_exists = false;
    for (CDT::Vertex_handle v : cdt.all_vertex_handles())
        if (v->point() == new_point) {
            point_exists = true;
            break;
        }
    if (point_exists)
        return true;

    if (!inside_boundaries(input1, new_point, points))
        return true;

    bool no_flip = false;
    CDT cdt_copy(cdt);
    cdt_copy.insert(new_point);
    if (count_obtuse_faces(cdt_copy) > count_obtuse_faces(cdt)) {
        cdt_copy = cdt;
        cdt_copy.insert_no_flip(new_point);
        if (count_obtuse_faces(cdt_copy) > count_obtuse_faces(cdt))
            return true;
        no_flip = true;
    }
    if (no_flip)
        cdt.insert_no_flip(new_point);
    else
        cdt.insert(new_point);
    return false;
}

int add_steiner_points_center(CDT& cdt, const int max_limit, Input1& input1, vector<Point>& points) {
    int steiner_points = 0;
    for (Face_handle_iterator f = cdt.finite_faces_begin(); f != cdt.finite_faces_end(); f++)
    {
        if (is_obtuse(cdt, f)) {
            if (steiner_method_center(cdt, input1, points, f))
                continue;

            steiner_points++;
            convergence.push_back(make_pair(cdt.finite_vertex_handles().size() - points.size(), count_obtuse_faces(cdt)));
            if (steiner_points > max_limit)
                return steiner_points;
            f = cdt.finite_faces_begin();
            for (int i = 0; i < rand() % cdt.number_of_faces() - 4; i++)
                f++;
        }
    }
    return steiner_points;
}

bool steiner_method_centroid(CDT &cdt, Input1 &input1, vector<Point> &points, Face_handle_iterator f) {
    Triangle t = cdt.triangle(f);
    Point p1 = f->vertex(0)->point();
    Point p2 = f->vertex(1)->point();
    Point p3 = f->vertex(2)->point();

    Point new_point = CGAL::centroid(t);

    bool point_exists = false;
    for (CDT::Vertex_handle v : cdt.all_vertex_handles())
        if (v->point() == new_point) {
            point_exists = true;
            break;
        }
    if (point_exists)
        return true;

    if (!inside_boundaries(input1, new_point, points))
        return true;

    bool no_flip = false;
    CDT cdt_copy(cdt);
    cdt_copy.insert(new_point);
    if (count_obtuse_faces(cdt_copy) > count_obtuse_faces(cdt)) {
        cdt_copy = cdt;
        cdt_copy.insert_no_flip(new_point);
        if (count_obtuse_faces(cdt_copy) > count_obtuse_faces(cdt))
            return true;
        no_flip = true;
    }
    if (no_flip)
        cdt.insert_no_flip(new_point);
    else
        cdt.insert(new_point);
    return false;
}

int add_steiner_points_centroid(CDT& cdt, const int max_limit, Input1& input1, vector<Point>& points) {
    int steiner_points = 0;
    for (Face_handle_iterator f = cdt.finite_faces_begin(); f != cdt.finite_faces_end(); f++)
    {
        if (is_obtuse(cdt, f)) {
            if (steiner_method_centroid(cdt, input1, points, f))
                continue;

            steiner_points++;
            convergence.push_back(make_pair(cdt.finite_vertex_handles().size() - points.size(), count_obtuse_faces(cdt)));
            //cout << "Inserted (center): " << new_point << endl;
            //cout << "Average area per triangle: " << average_area(cdt) << endl;
            if (steiner_points > max_limit)
                return steiner_points;
            f = cdt.finite_faces_begin();
            for (int i = 0; i < rand() % cdt.number_of_faces() - 4; i++)
                f++;
        }
    }
    return steiner_points;
}

Point find_middle_of_largest_edge(CDT& cdt, Face_handle f)
{
    Triangle t = cdt.triangle(f);
    Point p1 = f->vertex(0)->point();
    Point p2 = f->vertex(1)->point();
    Point p3 = f->vertex(2)->point();

    double dist12 = to_double(CGAL::squared_distance(p1, p2));
    double dist13 = to_double(CGAL::squared_distance(p1, p3));
    double dist23 = to_double(CGAL::squared_distance(p2, p3));

    if (dist12 >= dist13 && dist12 >= dist23) {
        //edge1-2 is largest
        return CGAL::midpoint(p1, p2);
    }
    else if (dist13 >= dist12 && dist13 >= dist23) {
        //edge1-3 is largest
        return CGAL::midpoint(p1, p3);
    }
    else {
        //edge2-3 is largest
        return CGAL::midpoint(p2, p3);
    }
    cout << "Error case!" << endl;
    return CGAL::midpoint(p1, p2);
}

Point find_projection_of_largest_edge(CDT& cdt, Face_handle f)
{
    Triangle t = cdt.triangle(f);
    Point p1 = f->vertex(0)->point();
    Point p2 = f->vertex(1)->point();
    Point p3 = f->vertex(2)->point();

    double dist12 = to_double(CGAL::squared_distance(p1, p2));
    double dist13 = to_double(CGAL::squared_distance(p1, p3));
    double dist23 = to_double(CGAL::squared_distance(p2, p3));

    if (dist12 >= dist13 && dist12 >= dist23) {
        //edge1-2 is largest
        return Line_2(p1, p2).projection(p3);
    }
    else if (dist13 >= dist12 && dist13 >= dist23) {
        //edge1-3 is largest
        return Line_2(p1, p3).projection(p2);
    }
    else {
        //edge2-3 is largest
        return Line_2(p2, p3).projection(p1);
    }
    cout << "Error case!" << endl;
    return Line_2(p1, p2).projection(p3);
}

bool steiner_method_middle(CDT &cdt, Input1 &input1, vector<Point> &points, Face_handle_iterator f) {
    Triangle t = cdt.triangle(f);
    Point p1 = f->vertex(0)->point();
    Point p2 = f->vertex(1)->point();
    Point p3 = f->vertex(2)->point();

    //find circumcenter
    Point new_point = find_middle_of_largest_edge(cdt, f);
    //cout << "To insert: " << p1 << endl << p2 << endl << p3 << endl << new_point << endl;

    bool point_exists = false;
    for (CDT::Vertex_handle v : cdt.all_vertex_handles())
        if (v->point() == new_point) {
            point_exists = true;
            break;
        }
    if (point_exists)
        return true;

    if (!inside_boundaries(input1, new_point, points))
        return true;

    CDT cdt_copy(cdt);
    cdt_copy.insert_no_flip(new_point);
    if (count_obtuse_faces(cdt_copy) > count_obtuse_faces(cdt))
        return true;
    cdt.insert_no_flip(new_point);
    return false;
}

int add_steiner_points_middle(CDT& cdt, const int max_limit, Input1& input1, vector<Point>& points) {
    int steiner_points = 0;
    for (Face_handle_iterator f = cdt.finite_faces_begin(); f != cdt.finite_faces_end(); f++)
    {
        if (is_obtuse(cdt, f)) {
            if (steiner_method_middle(cdt, input1, points, f))
                continue;
            steiner_points++;
            convergence.push_back(make_pair(cdt.finite_vertex_handles().size() - points.size(), count_obtuse_faces(cdt)));
            if (steiner_points > max_limit)
                return steiner_points;
            f = cdt.finite_faces_begin();
            for (int i = 0; i < rand() % cdt.number_of_faces() - 4; i++)
                f++;

        }
    }
    return steiner_points;
}

bool steiner_method_projection(CDT &cdt, Input1 &input1, vector<Point> &points, Face_handle_iterator f) {
    Triangle t = cdt.triangle(f);
    Point p1 = f->vertex(0)->point();
    Point p2 = f->vertex(1)->point();
    Point p3 = f->vertex(2)->point();

    //find projection
    Point new_point = find_projection_of_largest_edge(cdt, f);

    bool point_exists = false;
    for (CDT::Vertex_handle v : cdt.all_vertex_handles())
        if (v->point() == new_point) {
            point_exists = true;
            break;
        }
    if (point_exists)
        return true;

    if (!inside_boundaries(input1, new_point, points))
        return true;

    CDT cdt_copy(cdt);
    cdt_copy.insert_no_flip(new_point);
    if (count_obtuse_faces(cdt_copy) > count_obtuse_faces(cdt))
        return true;
    cdt.insert_no_flip(new_point);
    return false;
}

int add_steiner_points_projection(CDT& cdt, const int max_limit, Input1& input1, vector<Point>& points) {
    int steiner_points = 0;
    for (Face_handle_iterator f = cdt.finite_faces_begin(); f != cdt.finite_faces_end(); f++)
    {
        if (is_obtuse(cdt, f)) {
            if (steiner_method_projection(cdt, input1, points, f))
                continue;
            steiner_points++;
            convergence.push_back(make_pair(cdt.finite_vertex_handles().size() - points.size(), count_obtuse_faces(cdt)));
            if (steiner_points > max_limit)
                return steiner_points;
            f = cdt.finite_faces_begin();
            for (int i = 0; i < rand() % cdt.number_of_faces() - 4; i++)
                f++;

        }
    }
    return steiner_points;
}

void flip_edges(CDT &cdt)
{
    for (Face_handle f : cdt.finite_face_handles())
    {
        //cout << "f:" << endl;
        //cout << f->vertex(0)->point() << " " << f->vertex(1)->point() << " " << f->vertex(2)->point() << endl;
        if (is_obtuse(cdt, f)) {
            for (int n = 0; n < 3; n++) {
                Face_handle ni = f->neighbor(n);
                //cout << ni->vertex(0)->point() << " " << ni->vertex(1)->point() << " " << ni->vertex(2)->point() << endl;
                if (!ni->is_valid())
                    continue;
                if (f->is_constrained(n))
                    continue;
                if (is_finite(cdt, ni) == false)
                    continue;
                Polygon_2 temp;
                for (int i = 0; i < 3; i++) {
                    const Point& p = f->vertex(i)->point();
                    temp.push_back(p);
                }
                for (int i = 0; i < 3; i++) {
                    if (f->has_vertex(ni->vertex(i)))
                        continue;
                    temp.push_back(ni->vertex(i)->point());
                }
                if (temp.is_convex() == false)
                    continue;
                if (temp.vertices().size() <= 3)
                    continue;
                //if (cdt.is_flipable(f, n))
                    cdt.flip(f, n);
            }
        }
    }
}

bool is_obtuse(const CDT &cdt, Face_handle fh)
{
    Triangle t = cdt.triangle(fh);
    Point p1 = t.vertex(0);
    Point p2 = t.vertex(1);
    Point p3 = t.vertex(2);
    return (angle(p1, p2, p3) == CGAL::OBTUSE ||
        angle(p2, p3, p1) == CGAL::OBTUSE ||
        angle(p3, p1, p2) == CGAL::OBTUSE);
}

bool is_finite(const CDT& cdt, Face_handle fh)
{
    for (Face_handle f : cdt.finite_face_handles())
        if (fh == f)
            return true;
    return false;
}

int count_obtuse_faces(const CDT& cdt)
{
    int count = 0;
    for (Face_handle f : cdt.finite_face_handles())
        if (is_obtuse(cdt, f))
            count++;
    return count;
}

double get_rational(const K::FT& coord) {
    const auto exact_coord = CGAL::exact(coord);
    const mpq_t* gmpq_ptr = reinterpret_cast<const mpq_t*>(&exact_coord);

    mpz_t num, den;
    mpz_init(num);
    mpz_init(den);

    mpq_get_num(num, *gmpq_ptr);
    mpq_get_den(den, *gmpq_ptr);

    double res = mpz_get_ui(num)+0.0 / mpz_get_ui(den);

    mpz_clear(num);
    mpz_clear(den);

    // double num1 = atof(exact_coord.get_num().get_str().c_str());
    // double den1 = atof(exact_coord.get_den().get_str().c_str());
    return res;
}

Point get_centroid(const Polygon_2& polygon)
{
    double x = 0, y = 0;
    for (const Point &p: polygon) {
        x += get_rational(p.x());
        y += get_rational(p.y());
    }
    Point centroid(x / polygon.size(), y / polygon.size());
    return centroid;
}

double average_area(const CDT& cdt)
{
    double average = 0;
    for (Face_handle f : cdt.finite_face_handles())
        average += to_double(cdt.triangle(f).area());
    return average / cdt.finite_face_handles().size();
}

void print_result(string filename, Input2& input2, CDT& cdt, bool extras) {
    Input1& input = input2.input1;
    cout << "Writing output json file: " << filename << endl;
    ofstream out(filename);
    if (out.bad()) {
        perror("fopen");
        exit(1);
    }

    out << "{" << endl;
    out << quoted("content_type") << ": " << quoted("CG_SHOP_2025_Solution") << "," << endl;
    out << quoted("instance_uid") << ": " << quoted(input.instance_uid) << "," << endl;
    out << quoted("steiner_points_x") << ": [";
    int counter = 0;
    for (const Point& p : cdt.points()) {
        if (counter >= input.num_points) {
            out << "\""; print_rational(p.x(), out); out << "\"";
            if (counter != cdt.points().size() - 1)
                out << ", ";
        }
        counter++;
    }
    out << "],\n";
    out << quoted("steiner_points_y") << ": [";
    counter = 0;
    for (const Point& p : cdt.points()) {
        if (counter >= input.num_points) {
            out << "\""; print_rational(p.y(), out); out << "\"";
            if (counter != cdt.points().size() - 1)
                out << ", ";
        }
        counter++;
    }
    out << "],\n";
    out << quoted("edges") << ": [" << endl;
    counter = 0;
    for (const Edge& e : cdt.finite_edges()) {
        Point &p1 = e.first->vertex(e.first->cw(e.second))->point();
        Point &p2 = e.first->vertex(e.first->ccw(e.second))->point();
        out << "\t[" << locate_index(cdt, p1) << ", " << locate_index(cdt, p2) << "]";
        if (counter != cdt.finite_edges().size() - 1)
            out << ",";
        out << endl;
        counter++;
    }

    out << "],\n";
    out << quoted("obtuse_count") << ": " << count_obtuse_faces(cdt) << "," << endl;
    out << quoted("method") << ": " << quoted(input2.method) << "," << endl;
    out << quoted("parameters") << ": {";
    out << quoted("alpha") << ": " << input2.parameters.alpha << ", ";
    out << quoted("beta") << ": " << input2.parameters.beta << ", ";
    out << quoted("xi") << ": " << input2.parameters.xi << ", ";
    out << quoted("psi") << ": " << input2.parameters.psi << ", ";
    out << quoted("lambda") << ": " << input2.parameters.lambda << ", ";
    out << quoted("L") << ": " << input2.parameters.L << ", ";
    out << quoted("kappa") << ": " << input2.parameters.kappa << "},\n";
    out << quoted("randomization") << ": " << input2.randomization;
    if (extras)
        out << ",\n" << quoted("seed") << ": " << input2.seed;
    out << "\n}" << endl;
    out.close();
    cout << "Completed!" << endl;
}

void print_rational(const K::FT& coord, ofstream& fout) {
    const auto exact_coord = CGAL::exact(coord);
    //fout << exact_coord.get_num() << "/" << exact_coord.get_den();
    const mpq_t* gmpq_ptr = reinterpret_cast<const mpq_t*>(&exact_coord);

    mpz_t num, den;
    mpz_init(num);
    mpz_init(den);

    mpq_get_num(num, *gmpq_ptr);
    mpq_get_den(den, *gmpq_ptr);

    fout << num << "/" << den;

    mpz_clear(num);
    mpz_clear(den);
}

int locate_index(CDT& cdt, Point &p) {
    int pos = 0;
    for (Vertex_handle v: cdt.finite_vertex_handles()) {
        if (v->point() == p)
            return pos;
        pos++;
    }
    return -1;
}

bool inside_boundaries(Input1& input1, Point &p, vector<Point> points) {
    Polygon_2 polygon;
    for (int i = 0; i < input1.region_boundary.size(); i++)
        polygon.push_back(points[input1.region_boundary[i]]);
    if (polygon.is_simple() == false)
        return false;
    return polygon.bounded_side(p) != CGAL::ON_UNBOUNDED_SIDE;
}

//--------------------------------------------------------------------------

int LS_minimization(CDT& cdt, Input2& input, vector<Point>& points) {
    cout << "LS minimization" << endl;
    bool criterion = false;
    int current_iteration = 0;
    int previous_obtuse_count = -1;
    int steiner_points = 0;
    while(count_obtuse_faces(cdt) > 0 && criterion == false) {
        //for each obtuse triangle t do
        for (Face_handle_iterator f = cdt.finite_faces_begin(); f != cdt.finite_faces_end(); f++)
        {
            if (is_obtuse(cdt, f)) {
                int minimum = INT_MAX;
                CDT cdt_merge(cdt);
                CDT cdt_center(cdt);
                CDT cdt_centroid(cdt);
                CDT cdt_middle(cdt);
                CDT cdt_projection(cdt);
                CDT *best = NULL; // initialize best_T <- current triangulation

                if (!steiner_method_merge(cdt_merge, input.input1, points, f))
                    if (count_obtuse_faces(cdt_merge) < minimum) {
                        best = &cdt_merge;
                        minimum = count_obtuse_faces(cdt_merge);
                    }
                if (!steiner_method_center(cdt_center, input.input1, points, f))
                    if (count_obtuse_faces(cdt_center) < minimum) {
                        best = &cdt_center;
                        minimum = count_obtuse_faces(cdt_center);
                    }
                if (!steiner_method_centroid(cdt_centroid, input.input1, points, f))
                    if (count_obtuse_faces(cdt_centroid) < minimum) {
                        best = &cdt_centroid;
                        minimum = count_obtuse_faces(cdt_centroid);
                    }
                if (!steiner_method_middle(cdt_middle, input.input1, points, f))
                    if (count_obtuse_faces(cdt_middle) < minimum) {
                        best = &cdt_middle;
                        minimum = count_obtuse_faces(cdt_middle);
                    }
                if (!steiner_method_projection(cdt_projection, input.input1, points, f))
                    if (count_obtuse_faces(cdt_projection) < minimum) {
                        best = &cdt_projection;
                        minimum = count_obtuse_faces(cdt_projection);
                    }
                if (minimum != count_obtuse_faces(cdt) && best != NULL) {
                    cdt = *best;
                    convergence.push_back(make_pair(cdt.finite_vertex_handles().size() - points.size(), count_obtuse_faces(cdt)));
                    f = cdt.finite_faces_begin();
                }
            }
        }
        current_iteration++;

        // check stopping criterion
        if (previous_obtuse_count == count_obtuse_faces(cdt) || current_iteration == input.parameters.L)
            criterion = true;
    }
    steiner_points = cdt.all_vertex_handles().size() - points.size();
    cout << "LS minimization completed after " << current_iteration << " iterations." << endl;
    cout << "Inserted " << steiner_points << " steiner points" << endl;
    cout << "Final obtuse faces: " << count_obtuse_faces(cdt) << endl;
    CGAL::draw(cdt);
    return steiner_points;
}

double energy(CDT& cdt, Input2& input) {
    int steiner_points = cdt.all_vertex_handles().size() - input.input1.num_points;
    return input.parameters.alpha * count_obtuse_faces(cdt) + input.parameters.beta * steiner_points;
}

int simulated_annealing(CDT& cdt, Input2& input, vector<Point>& points) {
    cout << "Simulated Annealing" << endl;
    int steiner_points = 0;
    double T = 1; // initial temperature
    double previous_energy = energy(cdt, input);
    int initial_obtuse_faces = count_obtuse_faces(cdt);
    while (T >= 0) {
        cout << "T = " << T << endl;
        int iteration = 0;
        //for each obtuse triangle t do:
        for (Face_handle_iterator f = cdt.finite_faces_begin(); f != cdt.finite_faces_end(); f++)
        {
            if (is_obtuse(cdt, f)) {
                if (iteration++ > initial_obtuse_faces * 2)
                    break;
                //select randomly a steiner function
                bool (*steiner_method)(CDT&, Input1&, vector<Point>&, Face_handle_iterator);
                int myrand = rand() % 5;
                if (myrand == 0)
                    steiner_method = steiner_method_merge;
                else if (myrand == 1)
                    steiner_method = steiner_method_center;
                else if (myrand == 2)
                    steiner_method = steiner_method_centroid;
                else if (myrand == 3)
                    steiner_method = steiner_method_middle;
                else if (myrand == 4)
                    steiner_method = steiner_method_projection;

                CDT cdt_sa(cdt);
                if (!steiner_method(cdt_sa, input.input1, points, f)) {
                    convergence.push_back(make_pair(cdt.finite_vertex_handles().size() - points.size(), count_obtuse_faces(cdt)));
                    //inserted point, re-triangulated, calculating DE
                    double DE = energy(cdt_sa, input) - previous_energy;
                    double R = rand() / (double) RAND_MAX;
                    if (DE < 0 || exp(-DE / T) >= R) {
                        cdt = cdt_sa;
                        previous_energy = energy(cdt, input);
                        f = cdt.finite_faces_begin();
                    }
                }
            }
        }

        T = T - 1.0 / input.parameters.L;
    }
    steiner_points = cdt.all_vertex_handles().size() - points.size();
    cout << "Simulated Annealing completed" << endl;
    cout << "Inserted " << steiner_points << " steiner points" << endl;
    cout << "Final obtuse faces: " << count_obtuse_faces(cdt) << endl;
    CGAL::draw(cdt);
    return steiner_points;
}

double find_distance_of_largest_edge(CDT& cdt, Face_handle f)
{
    Triangle t = cdt.triangle(f);
    Point p1 = f->vertex(0)->point();
    Point p2 = f->vertex(1)->point();
    Point p3 = f->vertex(2)->point();

    double dist12 = to_double(CGAL::squared_distance(p1, p2));
    double dist13 = to_double(CGAL::squared_distance(p1, p3));
    double dist23 = to_double(CGAL::squared_distance(p2, p3));

    if (dist12 >= dist13 && dist12 >= dist23) {
        //edge1-2 is largest
        return sqrt(to_double(CGAL::squared_distance(p1, p2)));
    }
    else if (dist13 >= dist12 && dist13 >= dist23) {
        //edge1-3 is largest
        return sqrt(to_double(CGAL::squared_distance(p1, p3)));
    }
    else {
        //edge2-3 is largest
        return sqrt(to_double(CGAL::squared_distance(p2, p3)));
    }
    cout << "Error case!" << endl;
    return to_double(CGAL::squared_distance(p1, p2));;
}

double find_height_of_largest_edge(CDT& cdt, Face_handle f)
{
    Triangle t = cdt.triangle(f);
    Point p1 = f->vertex(0)->point();
    Point p2 = f->vertex(1)->point();
    Point p3 = f->vertex(2)->point();

    double dist12 = to_double(CGAL::squared_distance(p1, p2));
    double dist13 = to_double(CGAL::squared_distance(p1, p3));
    double dist23 = to_double(CGAL::squared_distance(p2, p3));

    if (dist12 >= dist13 && dist12 >= dist23) {
        //edge1-2 is largest
        Point proj = Line_2(p1, p2).projection(p3);
        return sqrt(to_double(CGAL::squared_distance(proj, p3)));
    }
    else if (dist13 >= dist12 && dist13 >= dist23) {
        //edge1-3 is largest
        Point proj = Line_2(p1, p3).projection(p2);
        return sqrt(to_double(CGAL::squared_distance(proj, p2)));
    }
    else {
        //edge2-3 is largest
        Point proj = Line_2(p2, p3).projection(p1);
        return sqrt(to_double(CGAL::squared_distance(proj, p1)));
    }
    cout << "Error case!" << endl;
    Point proj = Line_2(p1, p2).projection(p3);
    return to_double(CGAL::squared_distance(proj, p3));
}

double get_r(CDT& cdt, Face_handle f) {
    double R = find_distance_of_largest_edge(cdt, f) / 2;
    double h = find_height_of_largest_edge(cdt, f);
    // cout << f->vertex(0)->point() << ", " << f->vertex(1)->point() << ", " << f->vertex(2)->point() << endl
    //     << "\t" << R << ", " << h << endl;
    return R / h;
}

double heuristic_information_merge(CDT& cdt, Face_handle f) {
    return 1;
}

double heuristic_information_circumcenter(CDT& cdt, Face_handle f) {
    double r = get_r(cdt, f);
    return r / (2 + r);
}

double heuristic_information_centroid(CDT& cdt, Face_handle f) {
    return heuristic_information_circumcenter(cdt, f);
}

double heuristic_information_middle(CDT& cdt, Face_handle f) {
    double r = get_r(cdt, f);
    if (r < 0 || r >= 1.5)
        return 0;
    return (3 - 2 * r) / 3;
}

double heuristic_information_projection(CDT& cdt, Face_handle f) {
    double r = get_r(cdt, f);
    if (r < 0) return 0;
    return (r - 1) / r;
}

Face_handle_iterator random_picker(CDT& cdt) {
    int rand_num = rand() % cdt.finite_face_handles().size() + 1;
    for (Face_handle_iterator f = cdt.finite_faces_begin(); ; f++) {
        if (f == cdt.finite_faces_end())
            f = cdt.finite_faces_begin();
        if (is_obtuse(cdt, f))
            if (--rand_num == 0)
                return f;
    }
}

void* ant_function(void* args) {
    AntArgs& ant_args = *(AntArgs*)args;
    ant_args.selector = 0;
    srand(ant_args.ant_id * 1156 + ant_args.cycle_id * 87984 + ant_args.ant_id);
    //srand(pthread_self());
    //Improve triangulation(k)
    //randomly select an obtuse triangle
    if (count_obtuse_faces(ant_args.cdt_copy) == 0)
        pthread_exit(NULL);
    Face_handle_iterator random_face = random_picker(ant_args.cdt_copy);

    //calculate probabilities
    double probabilities[5] = {0, 0, 0, 0, 0}; //merge, circumcenter, centroid, middle, projection
    probabilities[0] = pow(ant_args.ti[0], ant_args.input.parameters.xi) * pow(heuristic_information_merge(ant_args.cdt_copy, random_face), ant_args.input.parameters.psi);
    probabilities[1] = pow(ant_args.ti[1], ant_args.input.parameters.xi) * pow(heuristic_information_circumcenter(ant_args.cdt_copy, random_face), ant_args.input.parameters.psi);
    probabilities[2] = pow(ant_args.ti[2], ant_args.input.parameters.xi) * pow(heuristic_information_centroid(ant_args.cdt_copy, random_face), ant_args.input.parameters.psi);
    probabilities[3] = pow(ant_args.ti[3], ant_args.input.parameters.xi) * pow(heuristic_information_middle(ant_args.cdt_copy, random_face), ant_args.input.parameters.psi);
    probabilities[4] = pow(ant_args.ti[4], ant_args.input.parameters.xi) * pow(heuristic_information_projection(ant_args.cdt_copy, random_face), ant_args.input.parameters.psi);
    double sum = 0;
    for (int i = 0; i < 5; i++)
        sum += probabilities[i];
    for (int i = 0; i < 5; i++)
        probabilities[i] /= sum;

    double myrand = rand() / (double) RAND_MAX;
    ant_args.selector = 0;
    double accumulator = 0;
    while(myrand >= accumulator) {
        accumulator += probabilities[ant_args.selector];
        if (myrand < accumulator)
            break;
        ant_args.selector++;
    }


    if (ant_args.selector == 0) {        //merge
        steiner_method_merge(ant_args.cdt_copy, ant_args.input.input1, ant_args.points, random_face);
    }
    else if (ant_args.selector == 1) {   //circumcenter
        steiner_method_center(ant_args.cdt_copy, ant_args.input.input1, ant_args.points, random_face);
    }
    else if (ant_args.selector == 2) {   //centroid
        steiner_method_centroid(ant_args.cdt_copy, ant_args.input.input1, ant_args.points, random_face);
    }
    else if (ant_args.selector == 3) {   //middle
        steiner_method_middle(ant_args.cdt_copy, ant_args.input.input1, ant_args.points, random_face);
    }
    else if (ant_args.selector >= 4) {   //projection
        steiner_method_projection(ant_args.cdt_copy, ant_args.input.input1, ant_args.points, random_face);
    }
    //return NULL;
    pthread_exit(NULL);
}

int ant_colony(CDT& cdt, Input2& input, vector<Point>& points) {
    cout << "Ant Colony" << endl;
    double t[5] = {0.8, 0.8, 0.8, 0.8, 0.8}; //merge, circumcenter, centroid, middle, projection
    for (int cycle = 0; cycle < input.parameters.L; cycle++) {
        pthread_t tid[input.parameters.kappa];
        AntArgs **ant_args = new AntArgs*[input.parameters.kappa];
        for (int i = 0; i < input.parameters.kappa; i++) {
            ant_args[i] = new AntArgs(i, cdt, input, t, points, cycle);
            pthread_create(&tid[i], NULL, ant_function, ant_args[i]);
            //ant_function(ant_args[i]);
        }
        for (int i = 0; i < input.parameters.kappa; i++)
            pthread_join(tid[i], NULL);

        //EvaluateTriangulation and SaveBest
        double min_energy = INT_MAX;
        double Dts[5] = {0, 0, 0, 0, 0};
        int pos = -1;
        for (int i = 0; i < input.parameters.kappa; i++) {
            double current_energy = energy(ant_args[i]->cdt_copy, input);
            if (count_obtuse_faces(ant_args[i]->cdt_copy) == 0)
                current_energy = 0;
            if (current_energy < min_energy) {
                min_energy = current_energy;
                pos = i;
            }
            if (current_energy != 0)
                Dts[ant_args[i]->selector] += 1 / (1 + current_energy);
        }
        CDT& best_cdt = ant_args[pos]->cdt_copy;

        //UpdatePheromones
        for (int i = 0; i < 5; i++)
            t[i] = (1 - input.parameters.lambda) * t[i] + Dts[i];

        cdt = best_cdt;
        convergence.push_back(make_pair(cdt.finite_vertex_handles().size() - points.size(), count_obtuse_faces(cdt)));
        for (int i = 0; i < input.parameters.kappa; i++)
            delete ant_args[i];
        delete[] ant_args;
    }
    int steiner_points = cdt.all_vertex_handles().size() - points.size() - 1;
    cout << "Ant Colony completed" << endl;
    cout << "Inserted " << steiner_points << " steiner points" << endl;
    cout << "Final obtuse faces: " << count_obtuse_faces(cdt) << endl;
    CGAL::draw(cdt);
    return steiner_points;
}

void remove_all_constraints(CDT& cdt) {
    // cout << "Constraints: " << cdt.constrained_edges().size() << endl;
    for (Face_handle f : cdt.finite_face_handles()) {
        for (int i = 0; i < 3; i++)
            cdt.remove_constrained_edge(f, i);
    }
    // cout << "Constraints: " << cdt.constrained_edges().size() << endl;
}

vector<pair<int, int>> save_all_constraints(CDT& cdt) {
    vector<pair<int, int>> constraints;
    for (Edge e : cdt.constrained_edges()) {
        Point &p1 = e.first->vertex(e.first->cw(e.second))->point();
        Point &p2 = e.first->vertex(e.first->ccw(e.second))->point();
        constraints.push_back(make_pair(locate_index(cdt, p1), locate_index(cdt, p2)));
    }
    return constraints;
}

void load_all_constraints(CDT& cdt, vector<pair<int, int>>& constraints, const vector<Point>& points) {
    for (const auto &constraint: constraints)
        cdt.insert_constraint(points[constraint.first], points[constraint.second]);
}

vector<pair<int, int>> polygon_to_constraints(CDT& cdt, Polygon_2& polygon) {
    vector<pair<int, int>> constraints;
    Polygon_2 convexHullPolygon;
    const Polygon_2::Vertices& range = polygon.vertices();
    vector<K::Point_2> res;
    CGAL::convex_hull_2(range.begin(), range.end(), back_inserter(res));
    for(auto it = res.begin(); it!= res.end(); ++it)
        convexHullPolygon.push_back(*it);

    for (Polygon_2::Edge_const_iterator e = convexHullPolygon.edges_begin(); e != convexHullPolygon.edges_end(); e++) {
        Point p1 = e->start();
        Point p2 = e->end();
        constraints.push_back(make_pair(locate_index(cdt, p1), locate_index(cdt, p2)));
    }
    return constraints;
}

//------------ part 3 -------------

Category get_CDT_Category(CDT& cdt, Input2& input, vector<Point>& points) {
    //Check convex
    Polygon_2 polygon;
    for (int i = 0; i < input.input1.region_boundary.size(); i++) {
        polygon.push_back(points[input.input1.region_boundary[i]]);
    }
    if (polygon.is_convex() == false) {
        //D or E
        if (input.input1.num_constraints > 0)
            return E;
        if (input.input1.num_constraints == 0) {
            for (int i = 1; i < input.input1.region_boundary.size(); i++)
                if (points[input.input1.region_boundary[i]].x() != points[input.input1.region_boundary[i - 1]].x() &&
                    points[input.input1.region_boundary[i]].y() != points[input.input1.region_boundary[i - 1]].y())
                    return E;
            return D;
        }
        return E;
    }
    else {
        //A, B or C
        if (input.input1.num_constraints == 0)
            return A;
        if (input.input1.num_constraints > 0) {
            for (int i = 1; i < input.input1.num_constraints; i++) {
                if (input.input1.additional_constraints[i].y != input.input1.additional_constraints[i - 1].x)
                    return B;
            }
            if (input.input1.additional_constraints[input.input1.num_constraints - 1].y !=
                input.input1.additional_constraints[0].x)
                return B;
            return C;
        }
    }

    return E;
}

bool steiner_method_random(CDT &cdt, Input1 &input1, vector<Point> &points, Face_handle_iterator f) {
    Triangle t = cdt.triangle(f);
    Point p1 = f->vertex(0)->point();
    Point p2 = f->vertex(1)->point();
    Point p3 = f->vertex(2)->point();
    //cout << p1 << " " << p2 << " " << p3 << endl;

    Point centr = centroid(p1, p2, p3), new_point;
    int tries = 50;
    for (int i = 0; i < tries; i++) {
        int iteration = 0, max_iterations = 3;
        do {
            if (iteration++ > max_iterations) {
                new_point = centr;
                break;
            }
            //Choose new random point
            int dist12 = sqrt(to_double(CGAL::squared_distance(p1, p2)));
            int dist13 = sqrt(to_double(CGAL::squared_distance(p1, p3)));
            int dist23 = sqrt(to_double(CGAL::squared_distance(p2, p3)));
            int dist, x = rand() % 3;
            if (x == 0)
                dist = dist12;
            else if (x == 1)
                dist = dist13;
            else if (x == 2)
                dist = dist23;
            if (dist == 0)
                dist = 1;
            int multiplier1 = rand() % 2 == 0 ? 1 : -1;
            int multiplier2 = rand() % 2 == 0 ? 1 : -1;
            Point p(new_point.x() + multiplier1 * rand() % dist, new_point.y() + multiplier2 * rand() % dist);
            new_point = p;
        } while(t.has_on_bounded_side(new_point) == false && t.has_on_boundary(new_point) == false);
        //cout << new_point << endl;

        bool point_exists = false;
        for (CDT::Vertex_handle v : cdt.all_vertex_handles())
            if (v->point() == new_point) {
                point_exists = true;
                break;
            }
        if (point_exists)
            return true;

        if (!inside_boundaries(input1, new_point, points))
            return true;

        bool to_insert = false;
        bool no_flip = false;
        CDT cdt_copy(cdt);
        cdt_copy.insert(new_point);
        if (count_obtuse_faces(cdt_copy) > count_obtuse_faces(cdt)) {
            to_insert = false;
            cdt_copy = cdt;
            cdt_copy.insert_no_flip(new_point);
            if (count_obtuse_faces(cdt_copy) > count_obtuse_faces(cdt))
                to_insert = false;
            else
                to_insert = true;
            no_flip = true;
        }
        else
            to_insert = true;
        if (to_insert) {
            if (no_flip)
                cdt.insert_no_flip(new_point);
            else
                cdt.insert(new_point);
            return false;
        }
    }

    return true;
}

int add_steiner_points_random(CDT& cdt, const int max_limit, Input1& input1, vector<Point>& points) {
    int steiner_points = 0;
    if (count_obtuse_faces(cdt) == 0)
        return 0;
    Face_handle_iterator f = random_picker(cdt);
    if (!steiner_method_random(cdt, input1, points, f)) {
        convergence.push_back(make_pair(cdt.finite_vertex_handles().size() - points.size(), count_obtuse_faces(cdt)));
        return 1;
    }
    return steiner_points;
}

double convergence_rate() {
    if (convergence.size() == 1)
        return 1;

    double pn = 0.0;
    int max = 1;
    for (int i = 1; i < convergence.size(); i++) {
        // pairs of <steiner_points, obtuse_faces>
        if (convergence[i].second != 0 && convergence[i - 1].first != 0 && convergence[i - 1].second != 0) {
            if (log((double)convergence[i].first / convergence[i - 1].first) == 0)
                continue;
            double pni = log((double)convergence[i].second / convergence[i - 1].second) / log((double)convergence[i].first / convergence[i - 1].first);
            pn += pni;
            max = convergence[i].first;
            //cout << convergence[i].first << "     \t" << convergence[i].second << "    \t" << pni << endl;
        }
    }
    if (max > 1)
        max--;
    if (max == 0)
        return 1;
    return -(pn / max);
}