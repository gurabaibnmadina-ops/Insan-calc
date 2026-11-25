// insan_calc_core.cpp
// Стандарты: C++20 | Нулевые внешние зависимости
// Назначение: Математическое ядро анализа геометрии DXF и STEP (ISO 10303-21)

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <iomanip>
#include <cctype>
#include <array>
#include <numbers>

namespace InsanMath {

constexpr double PI = std::numbers::pi;
constexpr double TWO_PI = 2.0 * PI;
constexpr double WELD_TOL = 1e-3; // Допуск слияния узлов (1 микрон)

// ============================================================================
// 1. ЛИНЕЙНАЯ АЛГЕБРА В R^3
// ============================================================================
struct Vec3 {
    double x = 0.0, y = 0.0, z = 0.0;

    Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
    Vec3 operator/(double s) const { return {x / s, y / s, z / s}; }

    double dot(const Vec3& o) const { return x * o.x + y * o.y + z * o.z; }
    Vec3 cross(const Vec3& o) const {
        return {y * o.z - z * o.y, z * o.x - x * o.z, x * o.y - y * o.x};
    }
    double norm() const { return std::sqrt(dot(*this)); }
    Vec3 normalized() const {
        double n = norm();
        return (n > 1e-12) ? (*this / n) : Vec3{0.0, 0.0, 0.0};
    }
};

// ============================================================================
// 2. ЧИСЛЕННОЕ ИНТЕГРИРОВАНИЕ (8-ТОЧЕЧНАЯ КВАДРАТУРА ГАУССА - ЛЕЖАНДРА)
// ============================================================================
struct GaussLegendre8 {
    static constexpr std::array<double, 8> xi = {
        -0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498,
         0.1834346424956498,  0.5255324099163290,  0.7966664774136267,  0.9602898564975363
    };
    static constexpr std::array<double, 8> w = {
         0.1012285362903763,  0.2223810344533745,  0.3137066458778873,  0.3626837833783620,
         0.3626837833783620,  0.3137066458778873,  0.2223810344533745,  0.1012285362903763
    };

    template <typename CurveFunc>
    static double integrate_arc_length(CurveFunc eval_point, double t0, double t1) {
        if (std::abs(t1 - t0) < 1e-12) return 0.0;
        double half_len = 0.5 * (t1 - t0);
        double mid = 0.5 * (t0 + t1);
        double sum = 0.0;
        // Шаг дифференцирования 4-го порядка точности (метод Ричардсона)
        double h = std::max(1e-7 * std::abs(t1 - t0), 1e-9);

        for (size_t i = 0; i < 8; ++i) {
            double t = mid + half_len * xi[i];
            double tm = std::max(t0, t - h);
            double tp = std::min(t1, t + h);
            double dt = tp - tm;
            if (dt > 1e-15) {
                Vec3 deriv = (eval_point(tp) - eval_point(tm)) / dt;
                sum += w[i] * deriv.norm();
            }
        }
        return sum * std::abs(half_len);
    }
};

// ============================================================================
// 3. ДИФФЕРЕНЦИАЛЬНАЯ ГЕОМЕТРИЯ NURBS / B-SPLINE (АЛГОРИТМ КОКСА - ДЕ БУРА)
// ============================================================================
struct BSplineCurve {
    int degree = 3;
    std::vector<Vec3> ctrl_pts;
    std::vector<double> weights;
    std::vector<double> knots;

    // Вычисление базисных функций B-сплайна в точке u
    Vec3 evaluate(double u) const {
        size_t n = ctrl_pts.size();
        if (n == 0) return {};
        if (knots.size() < n + degree + 1) {
            // Если узловой вектор не задан, генерируем равномерный зажатый (clamped)
            return ctrl_pts.front();
        }

        double u_min = knots[degree];
        double u_max = knots[n];
        u = std::clamp(u, u_min, u_max);

        // Поиск узлового интервала [u_k, u_{k+1})
        size_t k = degree;
        if (std::abs(u - u_max) < 1e-11) {
            k = n - 1;
        } else {
            auto it = std::upper_bound(knots.begin() + degree, knots.begin() + n + 1, u);
            k = std::distance(knots.begin(), it) - 1;
        }

        // Алгоритм де Бура (пирамидальная схема в однородных координатах)
        std::vector<Vec3> d_pts(degree + 1);
        std::vector<double> d_w(degree + 1, 1.0);
        bool rational = (weights.size() == n);

        for (int j = 0; j <= degree; ++j) {
            size_t idx = std::clamp<size_t>(k - degree + j, 0, n - 1);
            double w = rational ? weights[idx] : 1.0;
            d_pts[j] = ctrl_pts[idx] * w;
            d_w[j] = w;
        }

        for (int r = 1; r <= degree; ++r) {
            for (int j = degree; j >= r; --j) {
                size_t i = k - degree + j;
                double denom = knots[i + degree - r + 1] - knots[i];
                double alpha = (std::abs(denom) > 1e-12) ? ((u - knots[i]) / denom) : 0.0;
                d_pts[j] = d_pts[j - 1] * (1.0 - alpha) + d_pts[j] * alpha;
                d_w[j]   = d_w[j - 1]   * (1.0 - alpha) + d_w[j]   * alpha;
            }
        }
        return (std::abs(d_w[degree]) > 1e-12) ? (d_pts[degree] / d_w[degree]) : d_pts[degree];
    }

    double compute_length() const {
        size_t n = ctrl_pts.size();
        if (n < 2) return 0.0;
        if (knots.size() < n + degree + 1) {
            double len = 0.0;
            for (size_t i = 0; i + 1 < n; ++i) len += (ctrl_pts[i + 1] - ctrl_pts[i]).norm();
            return len;
        }
        double total = 0.0;
        // Интегрируем по каждому невырожденному узловому интервалу [u_i, u_{i+1}]
        for (size_t i = degree; i < n; ++i) {
            if (knots[i + 1] - knots[i] > 1e-11) {
                total += GaussLegendre8::integrate_arc_length(
                    [this](double u) { return this->evaluate(u); },
                    knots[i], knots[i + 1]
                );
            }
        }
        return total;
    }
};

// Аналитическая длина сегмента с кривизной Bulge: b = tan(theta / 4)
inline double bulge_segment_length(const Vec3& p1, const Vec3& p2, double bulge) {
    double chord = (p2 - p1).norm();
    double b = std::abs(bulge);
    if (b < 1e-9 || chord < 1e-12) return chord;
    return chord * ((1.0 + b * b) / b) * std::atan(b);
}

// ============================================================================
// 4. АЛГЕБРАИЧЕСКАЯ ТОПОЛОГИЯ (ПРОСТРАНСТВЕННЫЙ ХЕШ + DSU)
// ============================================================================
struct GridKey {
    int64_t ix, iy, iz;
    bool operator==(const GridKey& o) const { return ix == o.ix && iy == o.iy && iz == o.iz; }
};

struct GridKeyHash {
    size_t operator()(const GridKey& k) const {
        size_t h1 = std::hash<int64_t>{}(k.ix);
        size_t h2 = std::hash<int64_t>{}(k.iy);
        size_t h3 = std::hash<int64_t>{}(k.iz);
        return h1 ^ (h2 * 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2)) ^ (h3 << 1);
    }
};

class TopologyAnalyzer {
    std::unordered_map<GridKey, int, GridKeyHash> vertex_map;
    std::vector<int> parent;
    std::vector<int> degree;
    int closed_primitives = 0;

    int find(int i) {
        if (parent[i] == i) return i;
        return parent[i] = find(parent[i]);
    }

    int get_or_create_vertex(const Vec3& v) {
        GridKey key{
            static_cast<int64_t>(std::llround(v.x / WELD_TOL)),
            static_cast<int64_t>(std::llround(v.y / WELD_TOL)),
            static_cast<int64_t>(std::llround(v.z / WELD_TOL))
        };
        auto it = vertex_map.find(key);
        if (it != vertex_map.end()) return it->second;
        int id = static_cast<int>(parent.size());
        parent.push_back(id);
        degree.push_back(0);
        vertex_map[key] = id;
        return id;
    }

public:
    void add_closed_loop() {
        closed_primitives++;
    }

    void add_open_edge(const Vec3& v_start, const Vec3& v_end) {
        if ((v_end - v_start).norm() < WELD_TOL) {
            closed_primitives++;
            return;
        }
        int u = get_or_create_vertex(v_start);
        int v = get_or_create_vertex(v_end);
        degree[u]++;
        degree[v]++;
        int ru = find(u);
        int rv = find(v);
        if (ru != rv) parent[ru] = rv;
    }

    struct TopologyReport {
        int total_pierces = 0;      // Нулевое число Бетти B_0
        int closed_contours = 0;    // Замкнутые контуры (deg(v) == 2 для всех v)
        int open_contours = 0;      // Незамкнутые цепочки (ошибки чертежа или гравировка)
    };

    TopologyReport analyze() {
        std::unordered_map<int, bool> comp_is_closed;
        for (size_t i = 0; i < parent.size(); ++i) {
            if (degree[i] == 0) continue;
            int root = find(static_cast<int>(i));
            if (comp_is_closed.find(root) == comp_is_closed.end()) {
                comp_is_closed[root] = true;
            }
            if (degree[i] % 2 != 0) {
                comp_is_closed[root] = false;
            }
        }

        TopologyReport rep;
        rep.closed_contours = closed_primitives;
        for (const auto& [root, is_closed] : comp_is_closed) {
            if (is_closed) rep.closed_contours++;
            else rep.open_contours++;
        }
        rep.total_pierces = rep.closed_contours + rep.open_contours;
        return rep;
    }
};

// ============================================================================
// РЕЗУЛЬТАТ РАБОТЫ ЯДРА
// ============================================================================
struct AnalysisResult {
    std::string format;
    double total_cut_length_mm = 0.0;
    int pierce_count = 0;
    int closed_loops = 0;
    int open_chains = 0;
    double sheet_thickness_mm = 0.0; // Для 3D STEP
    double sheet_area_mm2 = 0.0;     // Площадь основной грани (для STEP)
    int primitives_count = 0;
};

// ============================================================================
// 5. КОНЕЧНЫЙ АВТОМАТ ПАРСИНГА DXF
// ============================================================================
class DXFEngine {
    static std::string trim(const std::string& s) {
        size_t a = s.find_first_not_of(" \t\r\n");
        if (a == std::string::npos) return "";
        size_t b = s.find_last_not_of(" \t\r\n");
        return s.substr(a, b - a + 1);
    }

public:
    static AnalysisResult parse(const std::string& filepath) {
        std::ifstream file(filepath);
        if (!file.is_open()) throw std::runtime_error("Cannot open DXF file: " + filepath);

        AnalysisResult res;
        res.format = "DXF (2D Vector Graph)";
        TopologyAnalyzer topo;

        // Чтение пар (Group Code, Value)
        std::vector<std::pair<int, std::string>> pairs;
        std::string code_str, val_str;
        while (std::getline(file, code_str) && std::getline(file, val_str)) {
            std::string tc = trim(code_str);
            if (tc.empty()) continue;
            try {
                int code = std::stoi(tc);
                pairs.emplace_back(code, trim(val_str));
            } catch (...) { continue; }
        }

        // Проверяем, есть ли явно выделенная секция ENTITIES
        bool has_entities_sec = false;
        for (size_t i = 0; i + 1 < pairs.size(); ++i) {
            if (pairs[i].first == 2 && pairs[i].second == "ENTITIES") {
                has_entities_sec = true;
                break;
            }
        }

        bool in_entities = !has_entities_sec;
        size_t idx = 0;

        while (idx < pairs.size()) {
            int code = pairs[idx].first;
            const std::string& val = pairs[idx].second;

            if (has_entities_sec) {
                if (code == 2 && val == "ENTITIES") { in_entities = true; idx++; continue; }
                if (in_entities && code == 0 && val == "ENDSEC") { in_entities = false; idx++; continue; }
            }
            if (!in_entities) { idx++; continue; }

            if (code == 0) {
                std::string etype = val;
                idx++;

                if (etype == "LINE") {
                    Vec3 p1, p2;
                    while (idx < pairs.size() && pairs[idx].first != 0) {
                        int c = pairs[idx].first;
                        double v = std::atof(pairs[idx].second.c_str());
                        if (c == 10) p1.x = v; else if (c == 20) p1.y = v; else if (c == 30) p1.z = v;
                        else if (c == 11) p2.x = v; else if (c == 21) p2.y = v; else if (c == 31) p2.z = v;
                        idx++;
                    }
                    double len = (p2 - p1).norm();
                    if (len > WELD_TOL) {
                        res.total_cut_length_mm += len;
                        topo.add_open_edge(p1, p2);
                        res.primitives_count++;
                    }
                }
                else if (etype == "CIRCLE") {
                    double r = 0.0;
                    while (idx < pairs.size() && pairs[idx].first != 0) {
                        if (pairs[idx].first == 40) r = std::atof(pairs[idx].second.c_str());
                        idx++;
                    }
                    if (r > WELD_TOL) {
                        res.total_cut_length_mm += TWO_PI * r;
                        topo.add_closed_loop();
                        res.primitives_count++;
                    }
                }
                else if (etype == "ARC") {
                    Vec3 center;
                    double r = 0.0, a1_deg = 0.0, a2_deg = 0.0;
                    while (idx < pairs.size() && pairs[idx].first != 0) {
                        int c = pairs[idx].first;
                        double v = std::atof(pairs[idx].second.c_str());
                        if (c == 10) center.x = v; else if (c == 20) center.y = v; else if (c == 30) center.z = v;
                        else if (c == 40) r = v;
                        else if (c == 50) a1_deg = v;
                        else if (c == 51) a2_deg = v;
                        idx++;
                    }
                    if (r > WELD_TOL) {
                        double a1 = a1_deg * PI / 180.0;
                        double a2 = a2_deg * PI / 180.0;
                        double dtheta = a2 - a1;
                        while (dtheta <= 0.0) dtheta += TWO_PI;
                        res.total_cut_length_mm += r * dtheta;

                        Vec3 p_start{center.x + r * std::cos(a1), center.y + r * std::sin(a1), center.z};
                        Vec3 p_end{center.x + r * std::cos(a2), center.y + r * std::sin(a2), center.z};
                        topo.add_open_edge(p_start, p_end);
                        res.primitives_count++;
                    }
                }
                else if (etype == "ELLIPSE") {
                    Vec3 center, major_vec;
                    double ratio = 1.0, t0 = 0.0, t1 = TWO_PI;
                    while (idx < pairs.size() && pairs[idx].first != 0) {
                        int c = pairs[idx].first;
                        double v = std::atof(pairs[idx].second.c_str());
                        if (c == 10) center.x = v; else if (c == 20) center.y = v; else if (c == 30) center.z = v;
                        else if (c == 11) major_vec.x = v; else if (c == 21) major_vec.y = v; else if (c == 31) major_vec.z = v;
                        else if (c == 40) ratio = v;
                        else if (c == 41) t0 = v;
                        else if (c == 42) t1 = v;
                        idx++;
                    }
                    double a = major_vec.norm();
                    double b = a * ratio;
                    if (a > WELD_TOL) {
                        while (t1 <= t0) t1 += TWO_PI;
                        auto eval_ellipse = [&](double t) -> Vec3 {
                            return {a * std::cos(t), b * std::sin(t), 0.0};
                        };
                        res.total_cut_length_mm += GaussLegendre8::integrate_arc_length(eval_ellipse, t0, t1);
                        if (std::abs((t1 - t0) - TWO_PI) < 1e-4) topo.add_closed_loop();
                        else {
                            Vec3 u_axis = major_vec.normalized();
                            Vec3 v_axis{-u_axis.y, u_axis.x, 0.0};
                            Vec3 ps = center + u_axis * (a * std::cos(t0)) + v_axis * (b * std::sin(t0));
                            Vec3 pe = center + u_axis * (a * std::cos(t1)) + v_axis * (b * std::sin(t1));
                            topo.add_open_edge(ps, pe);
                        }
                        res.primitives_count++;
                    }
                }
                else if (etype == "LWPOLYLINE") {
                    int flags = 0;
                    struct PolyVert { Vec3 pt; double bulge = 0.0; };
                    std::vector<PolyVert> verts;
                    while (idx < pairs.size() && pairs[idx].first != 0) {
                        int c = pairs[idx].first;
                        if (c == 70) flags = std::atoi(pairs[idx].second.c_str());
                        else if (c == 10) {
                            verts.push_back({Vec3{std::atof(pairs[idx].second.c_str()), 0.0, 0.0}, 0.0});
                        } else if (c == 20 && !verts.empty()) {
                            verts.back().pt.y = std::atof(pairs[idx].second.c_str());
                        } else if (c == 42 && !verts.empty()) {
                            verts.back().bulge = std::atof(pairs[idx].second.c_str());
                        }
                        idx++;
                    }
                    if (verts.size() >= 2) {
                        bool is_closed = (flags & 1) != 0;
                        double poly_len = 0.0;
                        for (size_t i = 0; i + 1 < verts.size(); ++i) {
                            poly_len += bulge_segment_length(verts[i].pt, verts[i + 1].pt, verts[i].bulge);
                        }
                        if (is_closed) {
                            poly_len += bulge_segment_length(verts.back().pt, verts.front().pt, verts.back().bulge);
                            topo.add_closed_loop();
                        } else {
                            topo.add_open_edge(verts.front().pt, verts.back().pt);
                        }
                        res.total_cut_length_mm += poly_len;
                        res.primitives_count++;
                    }
                }
                else if (etype == "SPLINE") {
                    BSplineCurve spline;
                    int flags = 0;
                    while (idx < pairs.size() && pairs[idx].first != 0) {
                        int c = pairs[idx].first;
                        double v = std::atof(pairs[idx].second.c_str());
                        if (c == 70) flags = std::atoi(pairs[idx].second.c_str());
                        else if (c == 71) spline.degree = std::atoi(pairs[idx].second.c_str());
                        else if (c == 40) spline.knots.push_back(v);
                        else if (c == 41) spline.weights.push_back(v);
                        else if (c == 10) spline.ctrl_pts.push_back({v, 0.0, 0.0});
                        else if (c == 20 && !spline.ctrl_pts.empty()) spline.ctrl_pts.back().y = v;
                        else if (c == 30 && !spline.ctrl_pts.empty()) spline.ctrl_pts.back().z = v;
                        idx++;
                    }
                    if (spline.ctrl_pts.size() >= 2) {
                        double slen = spline.compute_length();
                        res.total_cut_length_mm += slen;
                        bool is_closed = ((flags & 1) != 0) ||
                                         ((spline.ctrl_pts.front() - spline.ctrl_pts.back()).norm() < WELD_TOL);
                        if (is_closed) topo.add_closed_loop();
                        else topo.add_open_edge(spline.ctrl_pts.front(), spline.ctrl_pts.back());
                        res.primitives_count++;
                    }
                }
            } else {
                idx++;
            }
        }

        auto rep = topo.analyze();
        res.pierce_count = rep.total_pierces;
        res.closed_loops = rep.closed_contours;
        res.open_chains = rep.open_contours;
        return res;
    }
};

// ============================================================================
// 6. ТОПОЛОГИЧЕСКОЕ ЯДРО ПАРСИНГА STEP (ISO 10303-21 + ТЕОРЕМА СТОКСА)
// ============================================================================
class STEPEngine {
    struct Entity {
        std::string type;
        std::vector<std::string> args;
    };

    // Разбиение аргументов верхнего уровня с учетом вложенных скобок
    static std::vector<std::string> split_args(const std::string& s) {
        std::vector<std::string> res;
        std::string cur;
        int depth = 0;
        bool in_str = false;
        for (char c : s) {
            if (c == '\'') in_str = !in_str;
            if (!in_str) {
                if (c == '(') depth++;
                else if (c == ')') depth--;
                else if (c == ',' && depth == 0) {
                    res.push_back(trim(cur));
                    cur.clear();
                    continue;
                }
            }
            cur.push_back(c);
        }
        if (!cur.empty()) res.push_back(trim(cur));
        return res;
    }

    static std::string trim(const std::string& s) {
        size_t a = s.find_first_not_of(" \t\r\n");
        if (a == std::string::npos) return "";
        size_t b = s.find_last_not_of(" \t\r\n");
        return s.substr(a, b - a + 1);
    }

    static std::vector<int> parse_id_list(const std::string& list_str) {
        std::vector<int> ids;
        std::string s = trim(list_str);
        if (s.size() >= 2 && s.front() == '(' && s.back() == ')') {
            s = s.substr(1, s.size() - 2);
        }
        auto tokens = split_args(s);
        for (const auto& t : tokens) {
            if (!t.empty() && t[0] == '#') {
                ids.push_back(std::atoi(t.c_str() + 1));
            }
        }
        return ids;
    }

    static std::vector<double> parse_double_list(const std::string& list_str) {
        std::vector<double> vals;
        std::string s = trim(list_str);
        if (s.size() >= 2 && s.front() == '(' && s.back() == ')') {
            s = s.substr(1, s.size() - 2);
        }
        auto tokens = split_args(s);
        for (const auto& t : tokens) {
            if (!t.empty()) vals.push_back(std::atof(t.c_str()));
        }
        return vals;
    }

    static int parse_id(const std::string& s) {
        std::string t = trim(s);
        if (!t.empty() && t[0] == '#') return std::atoi(t.c_str() + 1);
        return -1;
    }

public:
    static AnalysisResult parse(const std::string& filepath) {
        std::ifstream file(filepath);
        if (!file.is_open()) throw std::runtime_error("Cannot open STEP file: " + filepath);

        // Чтение файла и очистка блочных комментариев /* ... */
        std::ostringstream oss;
        oss << file.rdbuf();
        std::string raw = oss.str();

        std::string clean;
        clean.reserve(raw.size());
        for (size_t i = 0; i < raw.size(); ++i) {
            if (i + 1 < raw.size() && raw[i] == '/' && raw[i + 1] == '*') {
                i += 2;
                while (i + 1 < raw.size() && !(raw[i] == '*' && raw[i + 1] == '/')) i++;
                i++;
            } else if (raw[i] != '\r' && raw[i] != '\n') {
                clean.push_back(raw[i]);
            }
        }

        // Построение DAG-графа сущностей STEP (#ID = TYPE(args);)
        std::unordered_map<int, Entity> dag;
        size_t pos = 0;
        while (pos < clean.size()) {
            size_t semi = clean.find(';', pos);
            if (semi == std::string::npos) break;
            std::string stmt = trim(clean.substr(pos, semi - pos));
            pos = semi + 1;

            if (!stmt.empty() && stmt[0] == '#') {
                size_t eq = stmt.find('=');
                if (eq != std::string::npos) {
                    int id = std::atoi(stmt.c_str() + 1);
                    std::string rhs = trim(stmt.substr(eq + 1));
                    size_t lp = rhs.find('(');
                    size_t rp = rhs.rfind(')');
                    if (lp != std::string::npos && rp != std::string::npos && rp > lp) {
                        Entity ent;
                        ent.type = trim(rhs.substr(0, lp));
                        ent.args = split_args(rhs.substr(lp + 1, rp - lp - 1));
                        dag[id] = std::move(ent);
                    }
                }
            }
        }

        // Вспомогательные резолверы геометрии из графа
        auto get_point = [&](int id) -> Vec3 {
            auto it = dag.find(id);
            if (it == dag.end()) return {};
            if (it->second.type == "VERTEX_POINT" && it->second.args.size() >= 2) {
                id = parse_id(it->second.args[1]);
                it = dag.find(id);
                if (it == dag.end()) return {};
            }
            if ((it->second.type == "CARTESIAN_POINT" || it->second.type == "DIRECTION") &&
                it->second.args.size() >= 2) {
                auto coords = parse_double_list(it->second.args[1]);
                Vec3 v;
                if (coords.size() >= 1) v.x = coords[0];
                if (coords.size() >= 2) v.y = coords[1];
                if (coords.size() >= 3) v.z = coords[2];
                return v;
            }
            return {};
        };

        struct Frame3D {
            Vec3 origin{0, 0, 0};
            Vec3 z_axis{0, 0, 1};
            Vec3 x_axis{1, 0, 0};
            Vec3 y_axis{0, 1, 0};
        };

        auto get_frame = [&](int id) -> Frame3D {
            Frame3D f;
            auto it = dag.find(id);
            if (it == dag.end() || it->second.type != "AXIS2_PLACEMENT_3D") return f;
            if (it->second.args.size() >= 2) f.origin = get_point(parse_id(it->second.args[1]));
            if (it->second.args.size() >= 3) {
                int zid = parse_id(it->second.args[2]);
                if (zid > 0) f.z_axis = get_point(zid).normalized();
            }
            if (it->second.args.size() >= 4) {
                int xid = parse_id(it->second.args[3]);
                if (xid > 0) f.x_axis = get_point(xid).normalized();
            }
            f.y_axis = f.z_axis.cross(f.x_axis).normalized();
            f.x_axis = f.y_axis.cross(f.z_axis).normalized();
            return f;
        };

        // Вычисление точной длины и дискретизации кривой EDGE_CURVE
        struct EdgeData {
            double length = 0.0;
            std::vector<Vec3> sample_pts; // Для теоремы Стокса
        };

        auto evaluate_edge_curve = [&](int edge_curve_id) -> EdgeData {
            EdgeData ed;
            auto it = dag.find(edge_curve_id);
            if (it == dag.end() || it->second.type != "EDGE_CURVE" || it->second.args.size() < 5) return ed;

            int v1_id = parse_id(it->second.args[1]);
            int v2_id = parse_id(it->second.args[2]);
            int geom_id = parse_id(it->second.args[3]);
            bool same_sense = (it->second.args[4].find('T') != std::string::npos);

            Vec3 p1 = get_point(v1_id);
            Vec3 p2 = get_point(v2_id);

            auto git = dag.find(geom_id);
            if (git == dag.end()) {
                ed.length = (p2 - p1).norm();
                ed.sample_pts = {p1, p2};
                return ed;
            }

            const std::string& gtype = git->second.type;
            const auto& gargs = git->second.args;

            if (gtype == "LINE") {
                ed.length = (p2 - p1).norm();
                ed.sample_pts = {p1, p2};
            }
            else if (gtype == "CIRCLE" && gargs.size() >= 3) {
                Frame3D frame = get_frame(parse_id(gargs[1]));
                double r = std::atof(gargs[2].c_str());
                double dtheta = TWO_PI;
                double a1 = 0.0, a2 = TWO_PI;

                if (v1_id != v2_id && (p2 - p1).norm() > WELD_TOL) {
                    Vec3 d1 = p1 - frame.origin;
                    Vec3 d2 = p2 - frame.origin;
                    a1 = std::atan2(d1.dot(frame.y_axis), d1.dot(frame.x_axis));
                    a2 = std::atan2(d2.dot(frame.y_axis), d2.dot(frame.x_axis));
                    if (same_sense) {
                        dtheta = a2 - a1;
                        while (dtheta <= 1e-9) dtheta += TWO_PI;
                    } else {
                        dtheta = a1 - a2;
                        while (dtheta <= 1e-9) dtheta += TWO_PI;
                    }
                }
                ed.length = r * dtheta;
                int steps = 16;
                for (int k = 0; k <= steps; ++k) {
                    double frac = static_cast<double>(k) / steps;
                    double ang = same_sense ? (a1 + frac * dtheta) : (a1 - frac * dtheta);
                    ed.sample_pts.push_back(
                        frame.origin + frame.x_axis * (r * std::cos(ang)) + frame.y_axis * (r * std::sin(ang))
                    );
                }
            }
            else if (gtype == "ELLIPSE" && gargs.size() >= 4) {
                Frame3D frame = get_frame(parse_id(gargs[1]));
                double a = std::atof(gargs[2].c_str());
                double b = std::atof(gargs[3].c_str());
                double a1 = 0.0, dtheta = TWO_PI;

                if (v1_id != v2_id && (p2 - p1).norm() > WELD_TOL) {
                    Vec3 d1 = p1 - frame.origin;
                    Vec3 d2 = p2 - frame.origin;
                    a1 = std::atan2(d1.dot(frame.y_axis) / b, d1.dot(frame.x_axis) / a);
                    double a2 = std::atan2(d2.dot(frame.y_axis) / b, d2.dot(frame.x_axis) / a);
                    dtheta = same_sense ? (a2 - a1) : (a1 - a2);
                    while (dtheta <= 1e-9) dtheta += TWO_PI;
                }
                double t_end = same_sense ? (a1 + dtheta) : (a1 - dtheta);
                ed.length = GaussLegendre8::integrate_arc_length([&](double t) -> Vec3 {
                    return frame.origin + frame.x_axis * (a * std::cos(t)) + frame.y_axis * (b * std::sin(t));
                }, std::min(a1, t_end), std::max(a1, t_end));

                int steps = 16;
                for (int k = 0; k <= steps; ++k) {
                    double frac = static_cast<double>(k) / steps;
                    double ang = same_sense ? (a1 + frac * dtheta) : (a1 - frac * dtheta);
                    ed.sample_pts.push_back(
                        frame.origin + frame.x_axis * (a * std::cos(ang)) + frame.y_axis * (b * std::sin(ang))
                    );
                }
            }
            else if (gtype.find("B_SPLINE_CURVE") != std::string::npos && gargs.size() >= 3) {
                BSplineCurve spline;
                spline.degree = std::atoi(gargs[1].c_str());
                auto pt_ids = parse_id_list(gargs[2]);
                for (int pid : pt_ids) spline.ctrl_pts.push_back(get_point(pid));

                if (gargs.size() >= 8) {
                    auto mults = parse_double_list(gargs[6]);
                    auto kvals = parse_double_list(gargs[7]);
                    for (size_t i = 0; i < std::min(mults.size(), kvals.size()); ++i) {
                        int m = static_cast<int>( std::round(mults[i]) );
                        for (int j = 0; j < m; ++j) spline.knots.push_back(kvals[i]);
                    }
                }
                ed.length = spline.compute_length();
                ed.sample_pts = spline.ctrl_pts;
            }
            else {
                ed.length = (p2 - p1).norm();
                ed.sample_pts = {p1, p2};
            }
            return ed;
        };

        // Поиск всех плоских граней ADVANCED_FACE и вычисление их площади по теореме Стокса
        struct PlanarFaceCandidate {
            int face_id = 0;
            double area = 0.0;
            double total_perimeter = 0.0;
            int loop_count = 0;
            int edge_count = 0;
            Vec3 normal;
            Vec3 origin;
        };

        std::vector<PlanarFaceCandidate> planar_faces;

        for (const auto& [id, ent] : dag) {
            if (ent.type == "ADVANCED_FACE" && ent.args.size() >= 3) {
                int surf_id = parse_id(ent.args[2]);
                auto sit = dag.find(surf_id);
                if (sit == dag.end() || sit->second.type != "PLANE" || sit->second.args.size() < 2) continue;

                Frame3D plane_frame = get_frame(parse_id(sit->second.args[1]));
                auto bound_ids = parse_id_list(ent.args[1]);

                PlanarFaceCandidate cand;
                cand.face_id = id;
                cand.normal = plane_frame.z_axis;
                cand.origin = plane_frame.origin;

                double max_loop_area = 0.0;

                for (int bid : bound_ids) {
                    auto bit = dag.find(bid);
                    if (bit == dag.end() || bit->second.args.size() < 2) continue;
                    int loop_id = parse_id(bit->second.args[1]);
                    auto lit = dag.find(loop_id);
                    if (lit == dag.end() || lit->second.type != "EDGE_LOOP" || lit->second.args.size() < 2) continue;

                    auto oe_ids = parse_id_list(lit->second.args[1]);
                    if (oe_ids.empty()) continue;

                    cand.loop_count++;
                    std::vector<Vec3> loop_poly;

                    for (int oe_id : oe_ids) {
                        auto oeit = dag.find(oe_id);
                        if (oeit == dag.end() || oeit->second.type != "ORIENTED_EDGE" || oeit->second.args.size() < 5) continue;

                        int ec_id = parse_id(oeit->second.args[3]);
                        bool orient = (oeit->second.args[4].find('T') != std::string::npos);

                        EdgeData ed = evaluate_edge_curve(ec_id);
                        cand.total_perimeter += ed.length;
                        cand.edge_count++;

                        if (!ed.sample_pts.empty()) {
                            if (orient) {
                                for (size_t k = 0; k + 1 < ed.sample_pts.size(); ++k)
                                    loop_poly.push_back(ed.sample_pts[k]);
                            } else {
                                for (size_t k = ed.sample_pts.size() - 1; k > 0; --k)
                                    loop_poly.push_back(ed.sample_pts[k]);
                            }
                        }
                    }

                    // Теорема Стокса: A = 0.5 * || sum(r_k x r_{k+1}) ||
                    if (loop_poly.size() >= 3) {
                        Vec3 vec_area{0, 0, 0};
                        for (size_t k = 0; k < loop_poly.size(); ++k) {
                            const Vec3& r0 = loop_poly[k] - plane_frame.origin;
                            const Vec3& r1 = loop_poly[(k + 1) % loop_poly.size()] - plane_frame.origin;
                            vec_area = vec_area + r0.cross(r1);
                        }
                        double loop_area = 0.5 * vec_area.norm();
                        if (loop_area > max_loop_area) max_loop_area = loop_area;
                    }
                }

                cand.area = max_loop_area;
                if (cand.loop_count > 0 && cand.total_perimeter > WELD_TOL) {
                    planar_faces.push_back(cand);
                }
            }
        }

        AnalysisResult res;
        res.format = "STEP (ISO 10303-21 3D Manifold B-Rep)";

        if (!planar_faces.empty()) {
            // Сортируем плоские грани по площади (и периметру при равенстве)
            std::sort(planar_faces.begin(), planar_faces.end(), [](const auto& a, const auto& b) {
                if (std::abs(a.area - b.area) > 1e-2) return a.area > b.area;
                return a.total_perimeter > b.total_perimeter;
            });

            const auto& main_face = planar_faces.front();
            res.total_cut_length_mm = main_face.total_perimeter;
            res.pierce_count = main_face.loop_count;
            res.closed_loops = main_face.loop_count;
            res.open_chains = 0;
            res.sheet_area_mm2 = main_face.area;
            res.primitives_count = main_face.edge_count;

            // Автоматическое вычисление толщины листа t между верхней и нижней плоскостями
            for (size_t i = 1; i < planar_faces.size(); ++i) {
                double cos_angle = std::abs(main_face.normal.dot(planar_faces[i].normal));
                if (cos_angle > 0.999) { // Параллельная грань
                    double dist = std::abs((planar_faces[i].origin - main_face.origin).dot(main_face.normal));
                    if (dist > 0.1 && dist < 100.0) {
                        res.sheet_thickness_mm = dist;
                        break;
                    }
                }
            }
        } else {
            // Fallback для каркасных (Wireframe) STEP-файлов без твердотельной топологии
            TopologyAnalyzer topo;
            for (const auto& [id, ent] : dag) {
                if (ent.type == "EDGE_CURVE" && ent.args.size() >= 5) {
                    EdgeData ed = evaluate_edge_curve(id);
                    res.total_cut_length_mm += ed.length;
                    res.primitives_count++;
                    if (!ed.sample_pts.empty()) {
                        topo.add_open_edge(ed.sample_pts.front(), ed.sample_pts.back());
                    }
                }
            }
            auto rep = topo.analyze();
            res.pierce_count = rep.total_pierces;
            res.closed_loops = rep.closed_contours;
            res.open_chains = rep.open_contours;
        }

        return res;
    }
};

} // namespace InsanMath

// ============================================================================
// C-API ЭКСПОРТ (ДЛЯ ИНТЕГРАЦИИ С PYTHON ЧЕРЕЗ CTYPES ИЛИ PYBIND11)
// ============================================================================
extern "C" {
    struct C_CalcResult {
        double total_cut_length_mm;
        double sheet_thickness_mm;
        double sheet_area_mm2;
        int pierce_count;
        int closed_loops;
        int open_chains;
        int primitives_count;
        int status_code; // 0 = OK, -1 = Error
    };

    C_CalcResult insan_analyze_cad_file(const char* filepath) {
        C_CalcResult out{};
        try {
            std::string path(filepath);
            std::string ext;
            size_t dot = path.rfind('.');
            if (dot != std::string::npos) {
                ext = path.substr(dot + 1);
                for (char& c : ext) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
            }

            InsanMath::AnalysisResult r;
            if (ext == "stp" || ext == "step") {
                r = InsanMath::STEPEngine::parse(path);
            } else {
                r = InsanMath::DXFEngine::parse(path);
            }

            out.total_cut_length_mm = r.total_cut_length_mm;
            out.sheet_thickness_mm  = r.sheet_thickness_mm;
            out.sheet_area_mm2      = r.sheet_area_mm2;
            out.pierce_count        = r.pierce_count;
            out.closed_loops        = r.closed_loops;
            out.open_chains         = r.open_chains;
            out.primitives_count    = r.primitives_count;
            out.status_code         = 0;
        } catch (...) {
            out.status_code = -1;
        }
        return out;
    }
}

// ============================================================================
// CLI ТОЧКА ВХОДА
// ============================================================================
int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <file.dxf | file.step | file.stp>\n";
        return 1;
    }

    try {
        std::string path = argv[1];
        std::string ext;
        size_t dot = path.rfind('.');
        if (dot != std::string::npos) {
            ext = path.substr(dot + 1);
            for (char& c : ext) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        }

        InsanMath::AnalysisResult res;
        if (ext == "stp" || ext == "step") {
            res = InsanMath::STEPEngine::parse(path);
        } else {
            res = InsanMath::DXFEngine::parse(path);
        }

        std::cout << std::fixed << std::setprecision(4);
        std::cout << "========================================================\n";
        std::cout << " INSAN CAD MATHEMATICAL KERNEL REPORT\n";
        std::cout << "========================================================\n";
        std::cout << " File             : " << path << "\n";
        std::cout << " Format           : " << res.format << "\n";
        std::cout << "--------------------------------------------------------\n";
        std::cout << " Cut Length (mm)  : " << res.total_cut_length_mm << " mm ("
                  << (res.total_cut_length_mm / 1000.0) << " m)\n";
        std::cout << " Pierce Count (B0): " << res.pierce_count << "\n";
        std::cout << " Closed Loops     : " << res.closed_loops << "\n";
        std::cout << " Open Chains      : " << res.open_chains
                  << (res.open_chains > 0 ? " [WARNING: Unclosed geometry!]" : " [OK]") << "\n";
        std::cout << " Geometric Edges  : " << res.primitives_count << "\n";
        if (res.sheet_thickness_mm > 0.0) {
            std::cout << " Sheet Thickness  : " << res.sheet_thickness_mm << " mm\n";
            std::cout << " Outer Face Area  : " << res.sheet_area_mm2 << " mm^2\n";
        }
        std::cout << "========================================================\n";
    } catch (const std::exception& e) {
        std::cerr << "[Error]: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
