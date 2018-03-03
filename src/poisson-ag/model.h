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
        geom::polygon < > m1, m2;
    };

    struct mesh_data
    {
        util::ptr_t < geom_data > geometry;
        util::ptr_t < geom::mesh > mesh;
        util::ptr_t < std::vector < double > > data;
        plot::triangulation_drawable :: ptr_t triangulation_plot;
        plot::dirichlet_cell_drawable :: ptr_t dirichlet_cell_plot;
        plot::drawable::ptr_t system_plot;
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
        n = size_t(std::floor(s / dl));
        for (size_t i = 0; i < n; ++i)
            path.points.emplace_back(-w / 2, -s / 2 - i * dl);
        n = size_t(std::floor(3 * w / 4 / dl));
        for (size_t i = 0; i < n; ++i)
            path.points.emplace_back(-w / 2 + i * dl, -h / 2);
        n = size_t(std::floor(M_PI * r1 / dl));
        for (size_t i = 0; i < n; ++i)
            path.points.emplace_back(w / 4 + r1 * std::sin(M_PI / n * i),
                                     - r1 * std::cos(M_PI / n * i));
        n = size_t(std::floor(3 * w / 4 / dl));
        for (size_t i = 0; i < n; ++i)
            path.points.emplace_back(w / 4 - i * dl, h / 2);
        n = size_t(std::floor(s / dl));
        for (size_t i = 0; i < n; ++i)
            path.points.emplace_back(-w / 2, h / 2 - i * dl);
        n = size_t(std::floor(3 * w / 4 / dl));
        for (size_t i = 0; i < n; ++i)
            path.points.emplace_back(-w / 2 + i * dl, s / 2);
        n = size_t(std::floor(M_PI * r2 / dl));
        for (size_t i = 0; i < n; ++i)
            path.points.emplace_back(w / 4 + r2 * std::sin(M_PI / n * i),
                                     r2 * std::cos(M_PI / n * i));
        n = size_t(std::floor(3 * w / 4 / dl));
        for (size_t i = 0; i < n; ++i)
            path.points.emplace_back(w / 4 - i * dl, -s / 2);

        return path;
    }

    inline geom::polygon < > transform_polygon(const geom::polygon < > & in,
                                               const geom::point2d_t & origin,
                                               double scale,
                                               double theta)
    {
        return geom::move(geom::scale(geom::rotate(in, theta), scale), origin);
    }

    inline geom_data make_geom(const parameters & p)
    {
        auto base = make_magnet_shape(min(p.dx, p.dy) / max(p.w, p.h) * 3);
        return
        {
            transform_polygon(base, p.m1_origin, p.m1_scale, p.m1_theta),
            transform_polygon(base, p.m2_origin, p.m2_scale, p.m2_theta)
        };
    }

    inline static plot::painter_t make_system_painter(const parameters & params,
                                                      mesh_data m)
    {
        return [&, m] (CDC & dc, const plot::viewport & vp)
        {
            auto metal_brush  = plot::palette::brush(RGB(0, 0, 0));
            auto border_brush = plot::palette::brush(RGB(0, 0, 0));
            auto m1_brush = plot::palette::brush(RGB(155, 0, 0));
            auto m2_brush = plot::palette::brush(RGB(0, 155, 0));

            auto border_pen = plot::palette::pen(0xac00ac, 3);

            auto m1_painter = plot::custom_drawable(plot::circle_painter(3, m1_brush));
            auto m2_painter = plot::custom_drawable(plot::circle_painter(3, m2_brush));
            auto metal_painter = plot::custom_drawable(plot::circle_painter(3, metal_brush));

            dc.SelectObject(border_pen.get());

            for (size_t i = 0, j = 1; i < m.geometry->m1.points.size(); ++i, ++j)
            {
                if (j == m.geometry->m1.points.size()) j = 0;
                dc.MoveTo(vp.world_to_screen().xy(m.geometry->m1.points[i]));
                dc.LineTo(vp.world_to_screen().xy(m.geometry->m1.points[j]));
            }

            for (size_t i = 0, j = 1; i < m.geometry->m2.points.size(); ++i, ++j)
            {
                if (j == m.geometry->m2.points.size()) j = 0;
                dc.MoveTo(vp.world_to_screen().xy(m.geometry->m2.points[i]));
                dc.LineTo(vp.world_to_screen().xy(m.geometry->m2.points[j]));
            }

            for (geom::mesh::idx_t i = 0; i < m.mesh->vertices().size(); ++i)
            {
                if (m.mesh->flags_at(i) & material::magnet1)
                    m1_painter.draw_at(dc, vp, m.mesh->point_at(i));
                else if (m.mesh->flags_at(i) & material::magnet2)
                    m2_painter.draw_at(dc, vp, m.mesh->point_at(i));
                else
                    metal_painter.draw_at(dc, vp, m.mesh->point_at(i));
            }
        };
    }

    inline mesh_data make_system_data(const parameters & p)
    {
        mesh_data md;

        std::vector < geom::point2d_t > super =
        {
            { -p.w, -p.h }, { p.w, -p.h }, { p.w, p.h }, { -p.w, p.h },
        };

        md.geometry = util::create < geom_data > (make_geom(p));
        auto & g = *md.geometry;

        md.mesh = util::create < geom::mesh > (false, false);

        md.mesh->init(super);

        md.mesh->add(g.m1.points, material::magnet1 | material::bound);
        md.mesh->add(g.m2.points, material::magnet2 | material::bound);

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
            if (geom::status::is(g.m1.contains(p0), geom::status::polygon::contains_point) ||
                geom::status::is(g.m2.contains(p0), geom::status::polygon::contains_point))
                continue;
            md.mesh->add(p0);
        }

        md.mesh->finish_mesh();

        md.data = util::create < std::vector < double > > (md.mesh->vertices().size());

        md.triangulation_plot = plot::triangulation_drawable::create(
            plot::make_data_source(md.mesh), nullptr);
        md.dirichlet_cell_plot = plot::dirichlet_cell_drawable::create(
            plot::make_data_source(md.mesh), nullptr, plot::palette::pen(RGB(155, 0, 0)));
        md.system_plot = plot::custom_drawable::create(
            make_system_painter(p, md));

        return md;
    }

    inline model_data make_model_data()
    {
        model_data md;
        md.config = make_plot_config();
        md.params = util::create < parameters > (make_default_parameters());
        md.field_line_data = make_plot_data();
        md.isoline_data = make_plot_data();
        md.system_data = make_system_data(*md.params);
        adjust(*md.params, *md.config.world);
        return md;
    }
}
