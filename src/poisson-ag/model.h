#pragma once

#include <afxwin.h>

#include <util/common/geom/geom.h>
#include <util/common/plot/plot.h>
#include <util/common/plot/triangulation_drawable.h>
#include <util/common/plot/dirichlet_cell_drawable.h>

namespace model
{

    using points_t = std::vector < geom::point2d_t > ;

    struct parameters
    {
        // system params
        double w, h;

        // geometry params
        double m1_scale, m1_theta, m2_scale, m2_theta;
        geom::point2d_t m1_origin, m2_origin;

        // other params
        double dx, dy, dxn, dyn;

        // material params
        double m1_qn, m1_qs, m2_qn, m2_qs;
    };

    inline static parameters make_default_parameters()
    {
        return
        {
            // system params
            100, 100,

            // geometry params
            30, 180, 30, 0,
            { -25, 0 }, { 25, 0 },

            // other params
            5, 5, 3, 3,

            // material params
            1, 1, -1, -1
        };
    }

    using material_t = geom::mesh::flags_t;

    namespace material
    {
        const material_t bound    = 0x1 << (0 + 10);
        const material_t ext      = 0x1 << (1 + 10);
        const material_t metal    = 0x1 << (2 + 10);
        const material_t dielectr = 0x1 << (3 + 10);
        const material_t magnet1  = 0x1 << (4 + 10);
        const material_t magnet2  = 0x1 << (5 + 10);
        const material_t north    = 0x1 << (6 + 10);
        const material_t south    = 0x1 << (7 + 10);
    };

    struct geom_data
    {
        geom::polygon < > m1_n, m1_s, m2_n, m2_s;
        geom::polygon < > m1_bb, m2_bb;
    };

    struct mesh_data
    {
        util::ptr_t < geom_data > geometry;
        util::ptr_t < geom::mesh > mesh;
        util::ptr_t < std::vector < double > > data;
        plot::triangulation_drawable :: ptr_t triangulation_plot;
        plot::dirichlet_cell_drawable :: ptr_t dirichlet_cell_plot;
        plot::drawable::ptr_t system_plot;
        plot::drawable::ptr_t point_plot;
    };

    struct plot_data
    {
        util::ptr_t < std::vector < points_t > > data;
        plot::multilist_drawable < points_t > :: ptr_t plot;
    };

    struct plot_config
    {
        plot::world_t::ptr_t world;
    };

    struct model_data
    {
        util::ptr_t < parameters > params;
        plot_config config;
        plot_data   isoline_data;
        plot_data   field_line_data;
        mesh_data   system_data;
    };

    inline static plot_data make_plot_data
    (
        plot::palette::pen_ptr pen = plot::palette::pen(0xffffff),
        plot::list_data_format data_format = plot::list_data_format::chain
    )
    {
        plot_data pd;
        pd.data = util::create < std::vector < points_t > > ();
        pd.plot = plot::multilist_drawable < points_t > :: create
        (
            plot::make_data_source(pd.data),
            nullptr, // no point painter
            pen
        );
        pd.plot->data_format = data_format;
        return pd;
    }

    inline static void adjust(const parameters & p,
                              plot::world_t & world)
    {
        world.xmin = - (world.xmax = p.w / 2);
        world.ymin = - (world.ymax = p.h / 2);
    }

    inline static plot::drawable::ptr_t make_root_drawable
    (
        const plot_config & p,
        std::vector < plot::drawable::ptr_t > layers
    )
    {
        using namespace plot;

        return viewporter::create(
            tick_drawable::create(
                layer_drawable::create(layers),
                const_n_tick_factory<axe::x>::create(
                    make_simple_tick_formatter(2, 5),
                    0,
                    5
                ),
                const_n_tick_factory<axe::y>::create(
                    make_simple_tick_formatter(2, 5),
                    0,
                    5
                ),
                palette::pen(RGB(80, 80, 80)),
                RGB(200, 200, 200)
            ),
            make_viewport_mapper(make_world_mapper(p.world))
        );
    }

    inline plot_config make_plot_config()
    {
        return { plot::world_t::create() };
    }

    inline geom::polygon < >
    make_magnet_shape(double dl)
    {
        const double w = 1, h = 0.5, s = h / 3,
                     r1 = h / 2, r2 = s / 2;

        geom::polygon < > path;

        size_t n;
        n = size_t(std::floor((w / 4 - r2) / dl));
        for (size_t i = 0; i < n; ++i)
            path.points.emplace_back(w / 4 + r2 + i * dl, 0);
        n = size_t(std::floor(M_PI / 2 * r1 / dl));
        for (size_t i = 0; i < n; ++i)
            path.points.emplace_back(w / 4 + r1 * std::cos(M_PI / 2 / n * i),
                                     + r1 * std::sin(M_PI / 2 / n * i));
        n = size_t(std::floor(3 * w / 4 / dl));
        for (size_t i = 0; i < n; ++i)
            path.points.emplace_back(w / 4 - i * dl, h / 2);
        n = size_t(std::floor(s / dl));
        for (size_t i = 0; i < n; ++i)
            path.points.emplace_back(-w / 2, h / 2 - i * dl);
        n = size_t(std::floor(3 * w / 4 / dl));
        for (size_t i = 0; i < n; ++i)
            path.points.emplace_back(-w / 2 + i * dl, s / 2);
        n = size_t(std::floor(M_PI / 2 * r2 / dl));
        for (size_t i = 0; i < n; ++i)
            path.points.emplace_back(w / 4 + r2 * std::sin(M_PI / 2 / n * i),
                                     r2 * std::cos(M_PI / 2 / n * i));

        return path;
    }

    inline geom::polygon < >
    make_magnet_bounding_box()
    {
        const double w = 1, h = 0.5;

        geom::polygon < > path;
        path.points.emplace_back(- w / 2, - h / 2);
        path.points.emplace_back(w / 2, - h / 2);
        path.points.emplace_back(w / 2, h / 2);
        path.points.emplace_back(- w / 2, h / 2);

        return path;
    }

    inline geom::polygon < > transform_polygon(const geom::polygon < > & in,
                                               const geom::point2d_t & origin,
                                               double scale_x, double scale_y,
                                               double theta)
    {
        auto p = in;
        for (size_t i = 0; i < p.points.size(); ++i)
        {
            auto p0 = geom::point2d_t(p.points[i].x * scale_x, p.points[i].y * scale_y);
            p.points[i] = p0.rotate(theta / 180 * M_PI) + origin;
        }
        return p;
    }

    inline geom_data make_geom(const parameters & p)
    {
        auto base = make_magnet_shape(min(p.dxn, p.dyn) / max(p.w, p.h));
        auto bb = make_magnet_bounding_box();
        return
        {
            transform_polygon(base, p.m1_origin, p.m1_scale, p.m1_scale, p.m1_theta),
            transform_polygon(base, p.m1_origin, p.m1_scale, -p.m1_scale, p.m1_theta),
            transform_polygon(base, p.m2_origin, p.m2_scale, p.m2_scale, p.m2_theta),
            transform_polygon(base, p.m2_origin, p.m2_scale, -p.m2_scale, p.m2_theta),
            transform_polygon(bb, p.m1_origin, 2 * p.m1_scale, 2 * p.m1_scale, p.m1_theta),
            transform_polygon(bb, p.m2_origin, 2 * p.m2_scale, 2 * p.m2_scale, p.m2_theta)
        };
    }

    inline static plot::painter_t make_system_painter(const parameters & params,
                                                      mesh_data m)
    {
        return [&, m] (CDC & dc, const plot::viewport & vp)
        {
            auto border_brush = plot::palette::brush(RGB(0, 0, 0));
            auto m1_brush = plot::palette::brush(RGB(155, 0, 0));
            auto m2_brush = plot::palette::brush(RGB(0, 155, 0));

            auto south_pen = plot::palette::pen(RGB(0, 0, 155), 5);
            auto north_pen = plot::palette::pen(RGB(155, 0, 0), 5);

            dc.SelectObject(south_pen.get());

            for (size_t i = 0, j = 1; i < m.geometry->m1_s.points.size(); ++i, ++j)
            {
                if (j == m.geometry->m1_s.points.size()) j = 0;
                dc.MoveTo(vp.world_to_screen().xy(m.geometry->m1_s.points[i]));
                dc.LineTo(vp.world_to_screen().xy(m.geometry->m1_s.points[j]));
            }

            for (size_t i = 0, j = 1; i < m.geometry->m2_s.points.size(); ++i, ++j)
            {
                if (j == m.geometry->m2_s.points.size()) j = 0;
                dc.MoveTo(vp.world_to_screen().xy(m.geometry->m2_s.points[i]));
                dc.LineTo(vp.world_to_screen().xy(m.geometry->m2_s.points[j]));
            }

            dc.SelectObject(north_pen.get());

            for (size_t i = 0, j = 1; i < m.geometry->m1_n.points.size(); ++i, ++j)
            {
                if (j == m.geometry->m1_n.points.size()) j = 0;
                dc.MoveTo(vp.world_to_screen().xy(m.geometry->m1_n.points[i]));
                dc.LineTo(vp.world_to_screen().xy(m.geometry->m1_n.points[j]));
            }

            for (size_t i = 0, j = 1; i < m.geometry->m2_n.points.size(); ++i, ++j)
            {
                if (j == m.geometry->m2_n.points.size()) j = 0;
                dc.MoveTo(vp.world_to_screen().xy(m.geometry->m2_n.points[i]));
                dc.LineTo(vp.world_to_screen().xy(m.geometry->m2_n.points[j]));
            }
        };
    }

    inline static plot::painter_t make_point_painter(const parameters & params,
                                                     mesh_data m)
    {
        return [&, m] (CDC & dc, const plot::viewport & vp)
        {
            auto metal_brush  = plot::palette::brush(RGB(0, 0, 0));

            auto point_painter = plot::custom_drawable(plot::circle_painter(3, metal_brush));

            for (geom::mesh::idx_t i = 0; i < m.mesh->vertices().size(); ++i)
            {
                point_painter.draw_at(dc, vp, m.mesh->point_at(i));
            }
        };
    }

    inline void update_system_data(const parameters & p, mesh_data & md)
    {
        std::vector < geom::point2d_t > super =
        {
            { -p.w, -p.h }, { p.w, -p.h }, { p.w, p.h }, { -p.w, p.h },
        };

        *md.geometry = make_geom(p);
        auto & g = *md.geometry;

        md.mesh->init(super);

        md.mesh->add(g.m1_s.points, material::magnet1 | material::south | material::bound);
        md.mesh->add(g.m1_n.points, material::magnet1 | material::north | material::bound);
        md.mesh->add(g.m2_s.points, material::magnet2 | material::south | material::bound);
        md.mesh->add(g.m2_n.points, material::magnet2 | material::north | material::bound);

        size_t n = size_t(std::floor(p.w / p.dx + 1));
        size_t m = size_t(std::floor(p.h / p.dy + 1));
        for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < m; ++j)
        {
            auto p0 = geom::make_point(-p.w / 2 + i * p.dx, -p.h / 2 + j * p.dy);
            if ((i == 0) || (i == (n - 1)) || (j == 0) || (j == (m - 1)))
            {
                md.mesh->add(p0, material::ext | material::bound);
                continue;
            }
            if (geom::status::is(g.m1_bb.contains(p0), geom::status::polygon::contains_point) ||
                geom::status::is(g.m2_bb.contains(p0), geom::status::polygon::contains_point))
                continue;
            p0.x += (rand() / (RAND_MAX + 1.) - 0.5) * p.dx / 5;
            p0.y += (rand() / (RAND_MAX + 1.) - 0.5) * p.dx / 5;
            md.mesh->add(p0);
        }

        double xmin = -p.w / 2,
               xmax = -p.w / 2 + (n - 1) * p.dx,
               ymin = -p.h / 2,
               ymax = -p.h / 2 + (m - 1) * p.dy;

        n = size_t(std::floor(p.w / p.dxn + 1));
        m = size_t(std::floor(p.h / p.dyn + 1));
        for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < m; ++j)
        {
            auto p0 = geom::make_point(-p.w / 2 + i * p.dxn, -p.h / 2 + j * p.dyn);
            if ((i == 0) || (i == (n - 1)) || (j == 0) || (j == (m - 1)))
                continue;
            if ((p0.x <= xmin) || (p0.x >= xmax) || (p0.y <= ymin) || (p0.y >= ymax))
                continue;
            if (geom::status::is_not(g.m1_bb.contains(p0), geom::status::polygon::contains_point) &&
                geom::status::is_not(g.m2_bb.contains(p0), geom::status::polygon::contains_point))
                continue;
            p0.x += (rand() / (RAND_MAX + 1.) - 0.5) * p.dxn / 5;
            p0.y += (rand() / (RAND_MAX + 1.) - 0.5) * p.dxn / 5;
            if (geom::status::is(g.m1_s.contains(p0), geom::status::polygon::contains_point) ||
                geom::status::is(g.m1_n.contains(p0), geom::status::polygon::contains_point) ||
                geom::status::is(g.m2_s.contains(p0), geom::status::polygon::contains_point) ||
                geom::status::is(g.m2_n.contains(p0), geom::status::polygon::contains_point))
                continue;
            md.mesh->add(p0);
        }

        md.mesh->finish_mesh();

        md.data = util::create < std::vector < double > > (md.mesh->vertices().size());
    }

    inline mesh_data make_system_data(const parameters & p)
    {
        mesh_data md;

        md.geometry = util::create < geom_data > ();
        md.mesh = util::create < geom::mesh > (false, false);

        update_system_data(p, md);

        md.triangulation_plot = plot::triangulation_drawable::create(
            plot::make_data_source(md.mesh), nullptr);
        md.dirichlet_cell_plot = plot::dirichlet_cell_drawable::create(
            plot::make_data_source(md.mesh), nullptr, plot::palette::pen(RGB(155, 0, 0)));
        md.system_plot = plot::custom_drawable::create(
            make_system_painter(p, md));
        md.point_plot = plot::custom_drawable::create(
            make_point_painter(p, md));

        return md;
    }

    inline model_data make_model_data(const parameters & p = make_default_parameters())
    {
        model_data md;
        md.config = make_plot_config();
        md.params = util::create < parameters > (p);
        md.field_line_data = make_plot_data();
        md.isoline_data = make_plot_data(plot::palette::pen(), plot::list_data_format::segment);
        md.system_data = make_system_data(*md.params);
        adjust(*md.params, *md.config.world);
        return md;
    }

    inline bool _get_val_intersection(double v1, double v2, double v, double & q)
    {
        double r = (v - v1) / (v2 - v1);
        if (isfinite(r) && (0 < r) && (r < 1))
        {
            q = r;
            return true;
        }
        return false;
    }

    inline void find_isolines(const geom::mesh & m,
                              const std::vector < double > & d,
                              double delta,
                              size_t max_lines,
                              std::vector < std::vector < geom::point2d_t > > & r)
    {
        r.clear();
        r.resize(max_lines * 2 + 1);
        for (geom::mesh::idx_t t = 0; t < m.triangles().size(); ++t)
        {
            auto & ti = m.triangles()[t];
            if (ti.flags & (geom::mesh::phantom | geom::mesh::superstruct)) continue;
            if (m.flags_at(ti.vertices[0]) &
                m.flags_at(ti.vertices[1]) &
                m.flags_at(ti.vertices[2]) &
                (material::magnet1 | material::magnet2)) continue;
            for (int l = - (int) max_lines; l <= (int) max_lines; ++l)
            {
                double val = l * delta;
                double q;
                size_t c = 0;
                if (_get_val_intersection(d[ti.vertices[0]], d[ti.vertices[1]], val, q))
                {
                    ++c;
                    r[(size_t) (l + max_lines)].push_back(
                        make_line_view(m.point_at(ti.vertices[0]),
                                       m.point_at(ti.vertices[1])).inner_point(q));
                }
                if (_get_val_intersection(d[ti.vertices[0]], d[ti.vertices[2]], val, q))
                {
                    ++c;
                    r[(size_t) (l + max_lines)].push_back(
                        make_line_view(m.point_at(ti.vertices[0]),
                                       m.point_at(ti.vertices[2])).inner_point(q));
                }
                if (_get_val_intersection(d[ti.vertices[1]], d[ti.vertices[2]], val, q))
                {
                    ++c;
                    r[(size_t) (l + max_lines)].push_back(
                        make_line_view(m.point_at(ti.vertices[1]),
                                       m.point_at(ti.vertices[2])).inner_point(q));
                }
                /* if occasionally added 1 or 3 points instead of 2;
                   may occur e.g. when the data contains !isfinite(data) */
                if ((c % 2) == 1)
                {
                    r[(size_t) (l + max_lines)].pop_back();
                }
            }
        }
    }

    /* *************************************************
       *********** Field Line Finder *******************
       ************************************************* */

    class field_line_finder
    {
    private:
        struct plane { double a, b, c; };
    private:
        util::ptr_t < geom::mesh > m;
        const std::vector < geom::mesh::idx_t > & hints;
        const std::vector < double > & d;
        std::vector < std::pair < bool, geom::point2d_t > > n;
        std::vector < bool > v;
    public:
        field_line_finder(util::ptr_t < geom::mesh > m,
                          const std::vector < double > & d,
                          const std::vector < geom::mesh::idx_t > & hints)
            : m(m)
            , d(d)
            , hints(hints)
        {
        }
    public:
        void find(std::vector < std::vector < geom::point2d_t > > & r);
    private:
        void _init();
        bool _normal(geom::mesh::idx_t t, geom::point2d_t & p) const;
        void _trace(geom::mesh::idx_t i, bool reverse, std::vector < geom::point2d_t > & r);
        bool _next(geom::mesh::idx_t i, geom::mesh::idx_t t, bool reverse, geom::point2d_t & p) const;
        bool _next(geom::mesh::idx_t & v0, geom::point2d_t & p, geom::mesh::idx_t & t, bool reverse) const;
        bool _next_triangle(geom::mesh::idx_t & v0, geom::mesh::idx_t & t0) const;
    };

    inline bool field_line_finder::_normal(
        geom::mesh::idx_t t, geom::point2d_t & p) const
    {
        auto & ti = m->triangles()[t];
        auto & p1 = m->point_at(ti.vertices[0]);
        auto & p2 = m->point_at(ti.vertices[1]);
        auto & p3 = m->point_at(ti.vertices[2]);

        double d0 = (p2.x - p3.x) * p1.y +
                    (p3.x - p1.x) * p2.y +
                    (p1.x - p2.x) * p3.y;

        double n1 = (p3.y - p2.y) * d[ti.vertices[0]] +
                    (p1.y - p3.y) * d[ti.vertices[1]] +
                    (p2.y - p1.y) * d[ti.vertices[2]];
        double n2 = (p2.x - p3.x) * d[ti.vertices[0]] +
                    (p3.x - p1.x) * d[ti.vertices[1]] +
                    (p1.x - p2.x) * d[ti.vertices[2]];

        p.x = n1 / d0;
        p.y = n2 / d0;

        auto n0 = math::norm(p);

        p.x /= n0;
        p.y /= n0;

        if (!isfinite(p.x) || !isfinite(p.y))
            return false;

        return true;
    }

    inline void field_line_finder::find(
        std::vector < std::vector < geom::point2d_t > > & r)
    {
        _init();
        r.resize(hints.size() * 2);
        for (size_t i = 0; i < hints.size(); ++i)
        {
            r[2 * i].clear();
            v.clear(); v.resize(m->triangles().size());
            r[2 * i].push_back(m->point_at(hints[i]));
            _trace(hints[i], false, r[2 * i]);
        }
        for (size_t i = 0; i < hints.size(); ++i)
        {
            r[2 * i + 1].clear();
            v.clear(); v.resize(m->triangles().size());
            r[2 * i + 1].push_back(m->point_at(hints[i]));
            _trace(hints[i], true, r[2 * i + 1]);
        }
    }

    inline void field_line_finder::_init()
    {
        n.resize(m->triangles().size());
        geom::point2d_t n0;
        for (geom::mesh::idx_t i = 0; i < m->triangles().size(); ++i)
        {
            if (m->triangles()[i].flags &
                (geom::mesh::phantom | geom::mesh::superstruct))
                continue;
            if (_normal(i, n0))
                n[i] = { true, n0 };
            else
                n[i].first = false;
        }
    }

    inline void field_line_finder::_trace(
        geom::mesh::idx_t i, bool reverse, std::vector < geom::point2d_t > & r)
    {
        geom::point2d_t p0;
        auto & nt = m->vertices()[i].neighbor_triangles;

        auto it = nt.begin();
        for (; it != nt.end(); ++it)
        {
            if (_next(i, *it, reverse, p0))
                break;
        }
        if (it == nt.end()) return;

        v[*it] = true;
        r.push_back(p0);

        geom::mesh::idx_t t = *it;
        geom::mesh::idx_t v0 = i;

        _next_triangle(v0, t);

        for (;;)
        {
            if (v[t])
                break;
            v[t] = true;
            if (!_next(v0, p0, t, reverse))
                break;
            r.push_back(p0);
            if ((v0 == SIZE_T_MAX) || (t == SIZE_T_MAX))
                break;
            if ((m->point_at(v0) == p0))
            {
                return _trace(v0, reverse, r);
            }
        }
    }

    inline bool field_line_finder::_next(
        geom::mesh::idx_t i, geom::mesh::idx_t t, bool reverse, geom::point2d_t & p) const
    {
        if (!std::get < 0 > (n[t])) return false;

        auto & ti = m->triangles()[t];

        geom::line l1, l2;

        l1.p1 = ((ti.vertices[0] != i) ?
                    m->point_at(ti.vertices[0]) :
                    m->point_at(ti.vertices[2]));
        l1.p2 = ((ti.vertices[1] != i) ?
                    m->point_at(ti.vertices[1]) :
                    m->point_at(ti.vertices[2]));
        l2.p1 = m->point_at(i);
        l2.p2 = m->point_at(i) + std::get < 1 > (n[t]);

        double q1, q2;
        auto s = l2.segment_intersection(l1, q1, q2);

        if (!geom::status::is_trusted(s, geom::status::line::intersects | geom::status::line::other_segment))
            return false;

        if ((q1 > 0) && !reverse || (q1 < 0) && reverse)
        {
            p = l2.inner_point(q1);
            return true;
        }

        return false;
    }


    inline bool field_line_finder::_next_triangle(
        geom::mesh::idx_t & v0, geom::mesh::idx_t & t0) const
    {
        auto & ti = m->triangles()[t0];

        auto v1 = (ti.vertices[0] == v0) ? ti.vertices[2] : ti.vertices[0];
        auto v2 = (ti.vertices[1] == v0) ? ti.vertices[2] : ti.vertices[1];

        v0 = SIZE_T_MAX;

        auto & nt = m->vertices()[v1].neighbor_triangles;

        for (auto it = nt.begin(); it != nt.end(); ++it)
        {
            if (*it == t0) continue;
            auto & ti0 = m->triangles()[*it];
            if ((ti0.vertices[0] == v1) && (ti0.vertices[1] == v2) ||
                (ti0.vertices[1] == v1) && (ti0.vertices[0] == v2))
                v0 = ti0.vertices[2];
            else if ((ti0.vertices[1] == v1) && (ti0.vertices[2] == v2) ||
                     (ti0.vertices[2] == v1) && (ti0.vertices[1] == v2))
                v0 = ti0.vertices[0];
            else if ((ti0.vertices[2] == v1) && (ti0.vertices[0] == v2) ||
                     (ti0.vertices[0] == v1) && (ti0.vertices[2] == v2))
                v0 = ti0.vertices[1];
            else continue;
            t0 = *it;
            break;
        }

        if (v0 == SIZE_T_MAX) t0 = SIZE_T_MAX;

        return v0 != SIZE_T_MAX;
    }

    inline bool field_line_finder::_next(
        geom::mesh::idx_t & v0, geom::point2d_t & p, geom::mesh::idx_t & t, bool reverse) const
    {
        if (!std::get < 0 > (n[t])) return false;

        auto & ti = m->triangles()[t];

        geom::line l1, l2;

        auto v1 = (ti.vertices[0] == v0) ? ti.vertices[2] : ti.vertices[0];
        auto v2 = (ti.vertices[1] == v0) ? ti.vertices[2] : ti.vertices[1];

        l1.p1 = m->point_at(v0);
        l1.p2 = m->point_at(v1);
        l2.p1 = p;
        l2.p2 = p + std::get < 1 > (n[t]);

        double q1, q2;

        auto s = l2.segment_intersection(l1, q1, q2);

        bool unsure = false;

        if (geom::status::is_trusted(s, geom::status::line::intersects))
        {
            auto t1 = geom::status::get(s, geom::status::line::other_segment);
            if ((t1 > 0) && ((q1 > 0) && !reverse || (q1 < 0) && reverse))
            {
                p = l1.inner_point(q2);
                v0 = v2;
                _next_triangle(v0, t);
                return true;
            }
            else if (t1 == 0)
            {
                unsure = true;
            }
        }

        l1.p2 = m->point_at(v2);
        l1.invalidate();

        s = l2.segment_intersection(l1, q1, q2);

        if (geom::status::is_trusted(s, geom::status::line::intersects))
        {
            auto t1 = geom::status::get(s, geom::status::line::other_segment);
            if ((t1 > 0) && ((q1 > 0) && !reverse || (q1 < 0) && reverse))
            {
                p = l1.inner_point(q2);
                v0 = v1;
                _next_triangle(v0, t);
                return true;
            }
            else if ((t1 == 0) && unsure && ((q1 > 0) && !reverse || (q1 < 0) && reverse))
            {
                p = m->point_at(v0);
                return true;
            }
        }

        return false;
    }

    /* *************************************************
       ******** Finite Element Method (Galerkin) *******
       ************************************************* */

    class finel_galerkin
    {
    private:
        struct plane { double a, b, c; };
        struct sparse_matrix
        {
            std::vector < std::vector < std::pair < size_t, double > > > matrix;
        };
    private:
        util::ptr_t < geom::mesh > m;
        /* here and below: a x = c */
        /* sparse matrix dim(a) = dim(x) * dim(x) */
        sparse_matrix a;
        /* dim(c) = dim(x) */
        std::vector < double > c;
        std::vector < double > x;
        /* mapping: variable -> mesh vertice */
        std::vector < geom::mesh::idx_t > vars;
        const size_t iters;
        const parameters & p;
    public:
        finel_galerkin(const parameters & p,
                       util::ptr_t < geom::mesh > m,
                       size_t iters = 100)
            : p(p)
            , m(m)
            , iters(iters)
        {
        }
    public:
        void next(std::vector < double > & r);
    private:
        void _init();
        void _next();
        plane _make_plane(geom::mesh::idx_t t,
                          geom::mesh::idx_t v) const;
        double _dot(geom::mesh::idx_t i,
                    geom::mesh::idx_t j) const;
        bool _is_var(geom::mesh::idx_t v) const;
        double _area(geom::mesh::idx_t t) const;
        double _charge_of(geom::mesh::idx_t i) const;
    };

    inline bool finel_galerkin::_is_var(geom::mesh::idx_t v) const
    {
        return (m->flags_at(v) &
                (material::ext | material::magnet1 | material::magnet2)) == 0;
    }

    inline double finel_galerkin::_charge_of(geom::mesh::idx_t i) const
    {
        if (m->flags_at(i) & material::north)
            if (m->flags_at(i) & material::magnet1)
                return p.m1_qn;
            else
                return p.m2_qn;
        else if (m->flags_at(i) & material::south)
            if (m->flags_at(i) & material::magnet1)
                return p.m1_qs;
            else
                return p.m2_qs;
        return 0;
    }

    inline double finel_galerkin::_area(geom::mesh::idx_t t) const
    {
        auto p1 = m->point_at(m->triangles()[t].vertices[1]) -
            m->point_at(m->triangles()[t].vertices[0]);
        auto p2 = m->point_at(m->triangles()[t].vertices[2]) -
            m->point_at(m->triangles()[t].vertices[0]);
        return std::abs(p1.x * p2.y - p1.y * p2.x) / 2;
    }

    inline finel_galerkin::plane finel_galerkin::_make_plane(
        geom::mesh::idx_t t, geom::mesh::idx_t v) const
    {
        auto & ti = m->triangles()[t];
        auto & p1 = m->point_at(ti.vertices[0]);
        auto & p2 = m->point_at(ti.vertices[1]);
        auto & p3 = m->point_at(ti.vertices[2]);
        double d = p1.x * (p2.y - p3.y) +
                   p2.x * (p3.y - p1.y) +
                   p3.x * (p1.y - p2.y);
        if (ti.vertices[0] == v)
        {
            return
            {
                (p2.y - p3.y) / d,
                (p3.x - p2.x) / d,
                (p2.x * p3.y - p3.x * p2.y) / d
            };
        }
        else if (ti.vertices[1] == v)
        {
            return
            {
                (p3.y - p1.y) / d,
                (p1.x - p3.x) / d,
                (p3.x * p1.y - p1.x * p3.y) / d
            };
        }
        else // if (ti.vertices[2] == v)
        {
            return
            {
                (p1.y - p2.y) / d,
                (p2.x - p1.x) / d,
                (p1.x * p2.y - p2.x * p1.y) / d
            };
        }
    }

    inline double finel_galerkin::_dot(
        geom::mesh::idx_t i, geom::mesh::idx_t j) const
    {
        double r = 0;
        auto & nt = m->vertices()[i].neighbor_triangles;
        for (auto it1 = nt.begin(); it1 != nt.end(); ++it1)
        {
            auto & ti = m->triangles()[*it1];
            if ((ti.vertices[0] != j) &&
                (ti.vertices[1] != j) &&
                (ti.vertices[2] != j))
                continue;
            auto p1 = _make_plane(*it1, i);
            for (auto it2 = nt.begin(); it2 != nt.end(); ++it2)
            {
                if (*it1 != *it2)
                    continue;
                auto p2 = _make_plane(*it1, j);
                r += - _area(*it1) * (
                    p1.a * p2.a +
                    p1.b * p2.b
                );
            }
        }
        return r;
    }

    inline void finel_galerkin::_init()
    {
        if (!vars.empty()) return;

        vars.resize(m->vertices().size());
        geom::mesh::idx_t var = 0;
        for (geom::mesh::idx_t i = 0; i < m->vertices().size(); ++i)
        {
            if (!_is_var(i)) continue;
            vars[var] = i;
            ++var;
        }
        vars.resize(var);

        x.resize(var);
	    x[0] = 0.5f;

        c.resize(var);
        a.matrix.resize(var);
        for (geom::mesh::idx_t j = 0, vj = 0; j < m->vertices().size(); ++j)
        {
            if (!_is_var(j)) continue;
            for (geom::mesh::idx_t i = 0, vi = 0; i < m->vertices().size(); ++i)
            {
                if (_is_var(i))
                {
                    auto d = _dot(i, j);
                    if (d != 0)
                        a.matrix[vj].emplace_back(vi, d);
                    ++vi;
                }
                else
                {
                    c[vj] += - _charge_of(i) * _dot(i, j);
                }
            }
            ++vj;
        }
    }

    inline void finel_galerkin::_next()
    {
        /* one iteration of Kaczmarz method
           without any accuracy checks */
        double s1, s2, t;
        for (size_t i = 0; i < a.matrix.size(); ++i)
        {
            s1 = s2 = 0;
		    for (size_t k = 0; k < a.matrix[i].size(); ++k)
		    {
                size_t j = std::get < 0 > (a.matrix[i][k]);
			    double v = std::get < 1 > (a.matrix[i][k]);
			    s1 += v * x[j];
			    s2 += v * v;
		    }
		    t = (c[i] - s1) / s2;
            for (size_t k = 0; k < a.matrix[i].size(); ++k)
		    {
                x[std::get < 0 > (a.matrix[i][k])] +=
                    std::get < 1 > (a.matrix[i][k]) * t;
		    }
        }
    }

    inline void finel_galerkin::next(std::vector < double > & r)
    {
        _init();
        for (size_t i = 0; i < iters; ++i) _next();
        r.resize(m->vertices().size());
        for (size_t j = 0, k = 0; k < m->vertices().size(); ++k)
        {
            if ((j < a.matrix.size()) && (k == vars[j]))
                r[k] = x[j++];
            else
                r[k] = _charge_of(k);
        }
    }
}
