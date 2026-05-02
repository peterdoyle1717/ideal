#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>

#include <array>
#include <cctype>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <fstream>
#include <vector>

namespace PMP = CGAL::Polygon_mesh_processing;

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point = Kernel::Point_3;
using Triangle = std::array<std::size_t, 3>;

[[noreturn]] static void die(const std::string& msg) {
    std::cerr << msg << '\n';
    std::exit(1);
}

static std::size_t parse_obj_index(const std::string& tok, std::size_t nverts) {
    const std::size_t slash = tok.find('/');
    const std::string head = tok.substr(0, slash);
    if (head.empty()) {
        throw std::runtime_error("bad OBJ face index");
    }
    std::size_t pos = 0;
    long long idx = 0;
    try {
        idx = std::stoll(head, &pos);
    } catch (...) {
        throw std::runtime_error("bad OBJ face index");
    }
    if (pos != head.size()) {
        throw std::runtime_error("bad OBJ face index");
    }
    if (idx > 0) {
        idx -= 1;
    } else if (idx < 0) {
        idx = static_cast<long long>(nverts) + idx;
    } else {
        throw std::runtime_error("OBJ indices are 1-based; got 0");
    }
    if (idx < 0 || static_cast<std::size_t>(idx) >= nverts) {
        throw std::runtime_error("OBJ face index out of range");
    }
    return static_cast<std::size_t>(idx);
}

static void read_obj(std::istream& in, std::vector<Point>& points, std::vector<Triangle>& triangles) {
    std::string line;
    while (std::getline(in, line)) {
        std::size_t k = 0;
        while (k < line.size() && std::isspace(static_cast<unsigned char>(line[k]))) {
            ++k;
        }
        if (k >= line.size() || line[k] == '#') {
            continue;
        }

        if (line[k] == 'v' && k + 1 < line.size() && std::isspace(static_cast<unsigned char>(line[k + 1]))) {
            std::istringstream iss(line.substr(k + 1));
            double x = 0.0, y = 0.0, z = 0.0;
            if (!(iss >> x >> y >> z)) {
                throw std::runtime_error("bad OBJ vertex line");
            }
            points.emplace_back(x, y, z);
            continue;
        }

        if (line[k] == 'f' && k + 1 < line.size() && std::isspace(static_cast<unsigned char>(line[k + 1]))) {
            std::istringstream iss(line.substr(k + 1));
            std::vector<std::string> toks;
            std::string tok;
            while (iss >> tok) {
                toks.push_back(tok);
            }
            if (toks.size() != 3) {
                throw std::runtime_error("expected triangular OBJ faces");
            }
            Triangle tri{
                parse_obj_index(toks[0], points.size()),
                parse_obj_index(toks[1], points.size()),
                parse_obj_index(toks[2], points.size())
            };
            triangles.push_back(tri);
            continue;
        }
    }

    if (points.empty()) {
        throw std::runtime_error("no vertices found");
    }
    if (triangles.empty()) {
        throw std::runtime_error("no faces found");
    }
}

static std::string stem(const std::string& path) {
    std::size_t slash = path.rfind('/');
    std::string base = (slash == std::string::npos) ? path : path.substr(slash + 1);
    std::size_t dot = base.rfind('.');
    if (dot != std::string::npos && base.substr(dot) == ".obj") base = base.substr(0, dot);
    return base;
}

static int check_stream(std::istream& in, const char* name) {
    try {
        std::vector<Point> points;
        std::vector<Triangle> triangles;
        read_obj(in, points, triangles);
        const bool ok = !PMP::does_triangle_soup_self_intersect<CGAL::Parallel_if_available_tag>(points, triangles);
        if (name)
            std::cout << (ok ? 1 : 0) << ' ' << name << '\n';
        else
            std::cout << (ok ? 1 : 0) << '\n';
        return ok ? 0 : 1;
    } catch (const std::exception& e) {
        std::cerr << "error";
        if (name) std::cerr << " (" << name << ")";
        std::cerr << ": " << e.what() << '\n';
        return 1;
    }
}

int main(int argc, char** argv) {
    if (argc < 2)
        return check_stream(std::cin, nullptr);

    int errors = 0;
    for (int i = 1; i < argc; i++) {
        std::ifstream in(argv[i]);
        if (!in) { std::cerr << "can't open: " << argv[i] << '\n'; errors++; continue; }
        errors += check_stream(in, stem(argv[i]).c_str());
    }
    return errors ? 1 : 0;
}
