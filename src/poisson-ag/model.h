#pragma once

#include <afxwin.h>

#include <map>

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
        double r;           /**< sphere radius           */
        geom::point2d_t dr; /**< charge displacement, [0,1]  */
        geom::point2d_t dr2; /**< charge displacement, [0,1]  */
        geom::point2d_t b;  /**< 2nd sphere displacement */

        // other params
        double s, ds;       /**< mean gap in the grid and
                                 the gap grow rate outside
                                 the sphere              */

        // material params
        double eps, q0, q02;
        BOOL has_q0, has_q02;
        BOOL has_1, has_2;
    };

    inline static parameters make_default_parameters()
    {
        return
        {
            // system params
            250, 250,

            // geometry params
            50, { 0, 0 }, { 0, 0 },
            { 120, 120 },

            // other params
            10, 1.3,

            // material params
            2, 1, -1,
            TRUE, FALSE,
            TRUE, TRUE
        };
    }

    using material_t = geom::mesh::flags_t;

    namespace material
    {
        const material_t bound    = 0x1 << (0 + 10);
        const material_t ext      = 0x1 << (1 + 10);
        const material_t metal    = 0x1 << (2 + 10);
        const material_t dielectr = 0x1 << (3 + 10);
        const material_t circle   = 0x1 << (4 + 10);
        const material_t charge_neighbor = 0x1 << (5 + 10);
        const material_t charge   = 0x1 << (6 + 10);
        const material_t charge2_neighbor = 0x1 << (7 + 10);
        const material_t charge2   = 0x1 << (8 + 10);
        const material_t circle_bound = 0x1 << (4 + 10);
    };

    struct geom_data
    {
        geom::polygon < > circle;
        geom::polygon < > circle2;
        std::vector < geom::mesh::idx_t > hints;
        std::map < geom::mesh::idx_t, std::pair < geom::mesh::idx_t, geom::mesh::idx_t > > bc_neighbors;
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
        world.xmin = - (world.xmax = p.w);
        world.ymin = - (world.ymax = p.h);
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
    make_circle_shape(double dl)
    {
        const double r = 1;

        geom::polygon < > path;

        size_t n = size_t(std::floor(2 * M_PI * r / dl));
        for (size_t i = 0; i < n; ++i)
            path.points.emplace_back(r * std::sin(2 * M_PI / n * i),
                                     r * std::cos(2 * M_PI / n * i));

        return path;
    }

    inline geom::polygon < > transform_polygon(const geom::polygon < > & in,
                                               const geom::point2d_t & origin,
                                               double scale)
    {
        auto p = in;
        for (size_t i = 0; i < p.points.size(); ++i)
        {
            p.points[i] = p.points[i] * scale + origin;
        }
        return p;
    }

    inline geom_data make_geom(const parameters & p)
    {
        auto base = make_circle_shape(p.s / p.r);
        return
        {
            transform_polygon(base, { 0, 0 }, p.r),
            transform_polygon(base, p.b, p.r)
        };
    }

    inline static plot::painter_t make_system_painter(const parameters & params,
                                                      mesh_data m)
    {
        return [&, m] (CDC & dc, const plot::viewport & vp)
        {
            auto border_brush = plot::palette::brush(RGB(0, 0, 0));

            auto circle_pen = plot::palette::pen(RGB(155, 155, 0), 5);

            dc.SelectObject(circle_pen.get());

            if (params.has_1)
            {
                for (size_t i = 0, j = 1; i < m.geometry->circle.points.size(); ++i, ++j)
                {
                    if (j == m.geometry->circle.points.size()) j = 0;
                    dc.MoveTo(vp.world_to_screen().xy(m.geometry->circle.points[i]));
                    dc.LineTo(vp.world_to_screen().xy(m.geometry->circle.points[j]));
                }
            }

            if (params.has_2)
            {
                for (size_t i = 0, j = 1; i < m.geometry->circle2.points.size(); ++i, ++j)
                {
                    if (j == m.geometry->circle2.points.size()) j = 0;
                    dc.MoveTo(vp.world_to_screen().xy(m.geometry->circle2.points[i]));
                    dc.LineTo(vp.world_to_screen().xy(m.geometry->circle2.points[j]));
                }
            }
        };
    }

    inline static plot::painter_t make_point_painter(const parameters & params,
                                                     mesh_data m)
    {
        return [&, m] (CDC & dc, const plot::viewport & vp)
        {
            auto metal_brush  = plot::palette::brush(0xacacac);

            auto point_painter = plot::custom_drawable(plot::circle_painter(2, metal_brush));

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
            { -2*p.w, -2*p.h }, { 2*p.w, -2*p.h }, { 2*p.w, 2*p.h }, { -2*p.w, 2*p.h },
        };

        *md.geometry = make_geom(p);
        auto & g = *md.geometry;

        md.mesh->init(super);

        double r0, s0;
        geom::point2d_t pdr = { p.r * p.dr.x, p.r * p.dr.y };
        double pdrn = math::norm(pdr);
        geom::point2d_t pdr2 = { p.r * p.dr2.x, p.r * p.dr2.y };
        double pdrn2 = math::norm(pdr2);

        if (p.has_1)
        {
            if (p.has_q0)
            {
                md.mesh->add(pdr, material::ext | material::bound | material::charge);

                s0 = p.s / 5;
                r0 = s0;
                {
                    auto cl = make_circle_shape(s0 / r0);
                    cl = transform_polygon(cl, pdr, r0);
                    md.mesh->add(cl.points.begin(), cl.points.end(), material::bound | material::charge_neighbor);
                }

                r0 = 2 * s0;
                while (r0 < min(p.s, p.r - pdrn))
                {
                    auto cl = make_circle_shape(s0 / r0);
                    cl = transform_polygon(cl, pdr, r0);
                    md.mesh->add(cl.points.begin(), cl.points.end());
                    r0 += s0;
                }
            }

            r0 = p.s;
            while (r0 < p.r - 1.5 * p.s)
            {
                auto cl = make_circle_shape(p.s / r0);
                cl = transform_polygon(cl, { 0, 0 }, r0);
                md.mesh->add(cl.points.begin(), cl.points.end());
                r0 += p.s;
            }

            r0 = p.r;
            {
                auto cl = make_circle_shape(p.s / r0);
                auto cm = transform_polygon(cl, { 0, 0 }, r0 - p.s);
                auto cp = transform_polygon(cl, { 0, 0 }, r0 + p.s);
                cl = transform_polygon(cl, { 0, 0 }, r0);
                auto idxm = md.mesh->add(cm.points.begin(), cm.points.end());
                auto idxl = md.mesh->add(cl.points.begin(), cl.points.end(), material::bound | material::circle_bound);
                auto idxp = md.mesh->add(cp.points.begin(), cp.points.end());
                for (size_t i = 0; i < idxl.size(); ++i)
                {
                    md.geometry->bc_neighbors[idxl[i]] = { idxm[i], idxp[i] };
                    g.hints.insert(g.hints.begin(), idxl.begin(), idxl.end());
                }
            }

            s0 = p.s;
            r0 = p.r + 2 * s0;
            while (r0 < max(p.w, p.h))
            {
                auto cl = make_circle_shape(s0 / r0);
                cl = transform_polygon(cl, { 0, 0 }, r0);
                md.mesh->add(cl.points.begin(), cl.points.end());
                r0 += s0;
                s0 += p.s * p.ds;
            }

            r0 = 1.5 * max(p.w, p.h);

            {
                auto cl = make_circle_shape(s0 / r0);
                cl = transform_polygon(cl, { 0, 0 }, r0);
                md.mesh->add(cl.points.begin(), cl.points.end(), material::ext | material::bound);
            }
        }

        /* other circle */

        if (p.has_2)
        {
            if (p.has_q02)
            {
                md.mesh->add(pdr2 + p.b, material::ext | material::bound | material::charge2);

                s0 = p.s / 5;
                r0 = s0;
                {
                    auto cl = make_circle_shape(s0 / r0);
                    cl = transform_polygon(cl, pdr2 + p.b, r0);
                    md.mesh->add(cl.points.begin(), cl.points.end(), material::bound | material::charge2_neighbor);
                }

                r0 = 2 * s0;
                while (r0 < min(p.s, p.r - pdrn2))
                {
                    auto cl = make_circle_shape(s0 / r0);
                    cl = transform_polygon(cl, pdr2 + p.b, r0);
                    md.mesh->add(cl.points.begin(), cl.points.end());
                    r0 += s0;
                }
            }

            md.mesh->add(p.b);

            r0 = p.s;
            while (r0 < p.r - 1.5 * p.s)
            {
                auto cl = make_circle_shape(p.s / r0);
                cl = transform_polygon(cl, p.b, r0);
                md.mesh->add(cl.points.begin(), cl.points.end());
                r0 += p.s;
            }

            r0 = p.r;
            {
                auto cl = make_circle_shape(p.s / r0);
                auto cm = transform_polygon(cl, p.b, r0 - p.s);
                auto cp = transform_polygon(cl, p.b, r0 + p.s);
                cl = transform_polygon(cl, p.b, r0);
                auto idxm = md.mesh->add(cm.points.begin(), cm.points.end());
                auto idxl = md.mesh->add(cl.points.begin(), cl.points.end(), material::bound | material::circle_bound);
                auto idxp = md.mesh->add(cp.points.begin(), cp.points.end());
                for (size_t i = 0; i < idxl.size(); ++i)
                {
                    md.geometry->bc_neighbors[idxl[i]] = { idxm[i], idxp[i] };
                }
            }

            if (!p.has_1)
            {
                s0 = p.s;
                r0 = p.r + 2 * s0;
                while (r0 < max(p.w, p.h))
                {
                    auto cl = make_circle_shape(s0 / r0);
                    cl = transform_polygon(cl, p.b, r0);
                    md.mesh->add(cl.points.begin(), cl.points.end());
                    r0 += s0;
                    s0 += p.s * p.ds;
                }

                r0 = 1.5 * max(p.w, p.h);

                {
                    auto cl = make_circle_shape(s0 / r0);
                    cl = transform_polygon(cl, p.b, r0);
                    md.mesh->add(cl.points.begin(), cl.points.end(), material::ext | material::bound);
                }
            }
        }

        // try to add q0 and q02 as regular points
        // may be helpful when the interaction energy
        // or something like it should be calculated
        md.mesh->add(pdr);
        md.mesh->add(pdr2 + p.b);

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
            plot::make_data_source(md.mesh), nullptr, plot::palette::pen(0xffffff));
        md.dirichlet_cell_plot = plot::dirichlet_cell_drawable::create(
            plot::make_data_source(md.mesh), nullptr, plot::palette::pen(RGB(255, 255, 0)));
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
        md.field_line_data = make_plot_data(plot::palette::pen(0x3cff3c, 2));
        md.isoline_data = make_plot_data(plot::palette::pen(0xffffff, 2), plot::list_data_format::segment);
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
                material::circle) continue;
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
        std::vector < geom::mesh::idx_t > vars_rev;
        const size_t iters;
        const parameters & p;
        const geom_data & gd;
    public:
        finel_galerkin(const parameters & p,
                       const geom_data & gd,
                       util::ptr_t < geom::mesh > m,
                       size_t iters = 100)
            : p(p)
            , gd(gd)
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
                (material::ext | material::bound)) == 0;
    }

    inline double finel_galerkin::_charge_of(geom::mesh::idx_t i) const
    {
        if (m->flags_at(i) & material::charge_neighbor)
            return p.q0 / p.eps / math::norm(geom::make_point(p.r * p.dr.x, p.r * p.dr.y) - m->point_at(i));
        if (m->flags_at(i) & material::charge2_neighbor)
            return p.q02 / p.eps / math::norm(geom::make_point(p.r * p.dr.x + p.b.x, p.r * p.dr.y + p.b.y) - m->point_at(i));
        if (m->flags_at(i) & material::circle_bound)
        {
            auto f1 = x[vars_rev[gd.bc_neighbors.at(i).first]],
                 f2 = x[vars_rev[gd.bc_neighbors.at(i).second]];
            return (p.eps * f1 + f2) / (p.eps + 1);
        }
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
        if (!vars.empty())
        {
            for (geom::mesh::idx_t j = 0, vj = 0; j < m->vertices().size(); ++j)
            {
                if (!_is_var(j)) continue;
                c[vj] = 0;
                for (geom::mesh::idx_t i = 0, vi = 0; i < m->vertices().size(); ++i)
                {
                    if (!_is_var(i))
                        c[vj] += - _charge_of(i) * _dot(i, j);
                }
                ++vj;
            }
            return;
        }

        vars.resize(m->vertices().size());
        vars_rev.resize(m->vertices().size(), -1);
        geom::mesh::idx_t var = 0;
        for (geom::mesh::idx_t i = 0; i < m->vertices().size(); ++i)
        {
            if (!_is_var(i)) continue;
            vars[var] = i;
            vars_rev[i] = var;
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
