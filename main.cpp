#include "includes/Tools.h"

using namespace std;

//Preselected parameter values
Parameters paramA(4.2, 1.3, 1, 1, 0.5, 10, 15);
Parameters paramB(4.2, 1.3, 1, 1, 0.5, 10, 15);
Parameters paramC(4.2, 1.3, 1, 1, 0.5, 10, 15);
Parameters paramD(3.0, 0.9, 1, 1, 0.5, 10, 20);
Parameters paramE(4.2, 1.3, 1, 1, 0.5, 10, 15);
extern vector<pair<int, int>> convergence;

int main(int argc, char** argv)
{
    string input_json, output_json;
    if (argc < 5 || argc > 8) {
        cout << "Usage: " << argv[0] << " -i <input_json> -o <output_json> -s <seed> -preselected_params" << endl;
        return 1;
    }
    input_json = argv[2];
    output_json = argv[4];
    Input2 input2;
    input2.seed = 0;
    bool preselected_params = false;
    for (int i = 1; i < argc; i++) {
        if (string(argv[i]) == "-preselected_params")
            preselected_params = true;
        else if (string(argv[i]) == "-s") {
            input2.seed = atoi(argv[i + 1]);
            i++;
        }
    }
    input2.parse(input_json);
    Input1& input1 = input2.input1;
    srand(input2.seed);

    // Initialize the Constrained Delaunay Triangulation (CDT)
    CDT cdt;

    // Define the points from the PSLG (x, y coordinates)
    vector<Point> points;
    for (int i = 0; i < input1.num_points; i++)
        points.push_back(Point(input1.points_x[i], input1.points_y[i]));

    // Insert points into the triangulation
    for (const Point &p: points)
        cdt.insert(p);

    // Define and add the constrained edges (from additional_constraints)
    vector<pair<int, int>> constraints;
    for (int i = 0; i < input1.num_constraints; i++)
        constraints.push_back({input1.additional_constraints[i].x, input1.additional_constraints[i].y});

    //Append also boundaries
    for (int i = 1; i < input1.region_boundary.size(); i++)
        constraints.push_back({input1.region_boundary[i], input1.region_boundary[i-1]});
    if (input1.region_boundary.size() > 1)
        constraints.push_back({input1.region_boundary[input1.region_boundary.size() - 1], input1.region_boundary[0]});

    // Insert constrained edges based on the provided indices
    for (const auto &constraint: constraints)
        cdt.insert_constraint(points[constraint.first], points[constraint.second]);

    int counter = 0;
    for (const Edge& e : cdt.finite_edges())
        if (cdt.is_constrained(e)) {
            counter++;
        }
    cout << "Number of constraints: " << counter << endl;
    cout << "Number of faces: " << cdt.number_of_faces() << endl;
    cout << "Finite edges: " << cdt.finite_edges().size() << endl;
    cout << "Initial obtuse faces: " << count_obtuse_faces(cdt) << endl;
    cout << "Average area per triangle: " << average_area(cdt) << endl;
    convergence.push_back(make_pair(0, count_obtuse_faces(cdt)));

    CGAL::draw(cdt);
    // vector<pair<int, int>> constraints_new = save_all_constraints(cdt);
    // remove_all_constraints(cdt);
    // CGAL::draw(cdt);
    // load_all_constraints(cdt, constraints_new, points);
    // CGAL::draw(cdt);

    Category category = get_CDT_Category(cdt, input2, points);
    if (category == A) {
        cout << "Category A" << endl;
        if (preselected_params) {
            input2.parameters = paramA;
            input2.method = "sa";
        }
    }
    else if (category == B) {
        cout << "Category B" << endl;
        if (preselected_params) {
            input2.parameters = paramB;
            //input2.parameters.kappa = points.size() / 4;
            input2.method = "ant";
        }
    }
    else if (category == C) {
        cout << "Category C" << endl;
        if (preselected_params) {
            input2.parameters = paramC;
            input2.method = "ant";
        }
    }
    else if (category == D) {
        cout << "Category D" << endl;
        if (preselected_params) {
            input2.parameters = paramD;
            input2.method = "sa";
            input2.delauney = true;
        }
    }
    else {
        cout << "Category E" << endl;
        if (preselected_params) {
            input2.parameters = paramE;
            input2.method = "sa";
        }
    }

    if (input2.method != "local" && input2.method != "sa" && input2.method != "ant") {
        cout << "method in " << argv[1] << " should be either local, sa or ant" << endl;
        return 1;
    }

    if (input2.delauney == false) {
        CDT copy_cdt(cdt);
        flip_edges(copy_cdt);
        if (count_obtuse_faces(copy_cdt) < count_obtuse_faces(cdt))
            cdt = copy_cdt;
        // Draw the triangulation using CGAL's draw function
        // CGAL::draw(cdt);
        cout << "After flip edge obtuse faces: " << count_obtuse_faces(cdt) << endl;

        for (int i = 0; i < points.size() / 10; i++) {
            cout << "Inserted " << add_steiner_points_merge(cdt, points.size() / 10, input1, points) << " steiner points using merge" << endl;
            cout << "Obtuse faces: " << count_obtuse_faces(cdt) << endl;
            cout << "Inserted " << add_steiner_points_center(cdt, points.size() / 10, input1, points) << " steiner points using center" << endl;
            cout << "Obtuse faces: " << count_obtuse_faces(cdt) << endl;
            cout << "Inserted " << add_steiner_points_centroid(cdt, points.size() / 10, input1, points) << " steiner points using centroid" << endl;
            cout << "Obtuse faces: " << count_obtuse_faces(cdt) << endl;
            cout << "Inserted " << add_steiner_points_projection(cdt, points.size() / 10, input1, points) << " steiner points using projection" << endl;
            cout << "Obtuse faces: " << count_obtuse_faces(cdt) << endl;
            cout << "Inserted " << add_steiner_points_middle(cdt, points.size() / 10, input1, points) << " steiner points using middle" << endl;
            cout << "Final obtuse faces: " << count_obtuse_faces(cdt) << endl;
            cout << "Obtuse faces: " << count_obtuse_faces(cdt) << endl;
        }
        CGAL::draw(cdt);
    }

    if (input2.method == "local")
        LS_minimization(cdt, input2, points);
    else if (input2.method == "sa")
        simulated_annealing(cdt, input2, points);
    else if (input2.method == "ant")
        ant_colony(cdt, input2, points);

    //Part 3
    int total_random = 0;
    int max_iterations = 0;
    if (points.size() > 50)
        max_iterations = points.size() / 5;
    else
        max_iterations = points.size();
    for (int i = 0; i < max_iterations; i++) {
        int inserted;
        if (inserted = add_steiner_points_random(cdt, max_iterations / 10, input1, points) > 0) {
            total_random += inserted;
        }
        //Print progress
        if (i > 0 && int(i * 100.0 / max_iterations) % 10 == 0)
            cout << i * 100.0 / max_iterations << "% completed..." << endl;
    }
    if (total_random > 0)
        input2.randomization = "true";
    cout << "Inserted " << total_random << " steiner points using random" << endl;
    cout << "Final obtuse faces: " << count_obtuse_faces(cdt) << endl;
    cout << "Convergence rate: " << convergence_rate() << endl;
    CGAL::draw(cdt);

    print_result(output_json, input2, cdt);

    return 0;
}

