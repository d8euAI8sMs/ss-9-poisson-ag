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
        double dx, dy;

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
            2, 2,

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
        auto base = make_magnet_shape(min(p.dx, p.dy) / max(p.w, p.h) * 3);
        return
        {
            transform_polygon(base, p.m1_origin, p.m1_scale, p.m1_scale, p.m1_theta),
            transform_polygon(base, p.m1_origin, p.m1_scale, -p.m1_scale, p.m1_theta),
            transform_polygon(base, p.m2_origin, p.m2_scale, p.m2_scale, p.m2_theta),
            transform_polygon(base, p.m2_origin, p.m2_scale, -p.m2_scale, p.m2_theta)
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

        size_t n = size_t(std::floor(p.w / p.dx));
        size_t m = size_t(std::floor(p.h / p.dy));
        for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < m; ++j)
        {
            auto p0 = geom::make_point(-p.w / 2 + i * p.dx, -p.h / 2 + j * p.dy);
            if ((i == 0) || (i == (n - 1)) || (j == 0) || (j == (m - 1)))
            {
                md.mesh->add(p0, material::ext | material::bound);
                continue;
            }
            p0.x += (rand() / (RAND_MAX + 1.) - 0.5) * p.dx / 5;
            p0.y += (rand() / (RAND_MAX + 1.) - 0.5) * p.dx / 5;
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
}
