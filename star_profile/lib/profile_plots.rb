# profiles.rb

class ProfilePlots

    include Math
    include Tioga
    include FigureConstants
    
    def background
        t.fill_color = OldLace
        t.fill_frame
    end
    
    def show_box_labels(title, xlabel=nil, ylabel=nil)
        if title != nil
            t.show_title(title); t.no_title
        end
        if xlabel != nil
            t.show_xlabel(xlabel); t.no_xlabel
        end
        if ylabel != nil
            t.show_ylabel(ylabel); t.no_ylabel
        end
    end
    
    def show_model_number(pos = 1.15, shift = 4)
        return if t.in_subplot or t.model_number < 0
        txt = sprintf('%i', t.model_number)
        t.show_text('text' => txt, 'side' => TOP, 'pos' => pos, 'shift' => shift, 'scale' => 0.8,
            'justification' => RIGHT_JUSTIFIED)
        t.show_text('text' => d.initial_Mass_String, 'side' => TOP, 'pos' => pos, 'shift' => shift-0.5, 'scale' => 0.5,
            'justification' => RIGHT_JUSTIFIED)
        t.show_text('text' => d.initial_Z_String, 'side' => TOP, 'pos' => pos, 'shift' => shift-1.9, 'scale' => 0.5,
            'justification' => RIGHT_JUSTIFIED)
        if d.star_Mass < 0.1
            mstr = '$M\!$=%0.4f'
        elsif d.star_Mass < 1
            mstr = '$M\!$=%0.3f'
        else
            mstr = '$M\!$=%0.2f'
        end
        t.show_text('text' => sprintf(mstr, d.star_Mass), 
            'side' => TOP, 'pos' => pos, 'shift' => shift-3.4, 'scale' => 0.5,
            'justification' => RIGHT_JUSTIFIED)
    end
    
    def get_star_age_title(fmt)
        return sprintf(fmt, d.star_Age_String)
    end
    
    def show_mass_locs(xs, ymin)
        return if xs == d.sx_M
        y1 = ymin + t.default_text_height_dy * 1.5
        t.save_legend_separator(1.5)
        x = xs[d.sx_M.where_closest(d.star_Mass * 0.5)]
        t.line_type = Line_Type_Solid
        t.line_width = 2
        t.stroke_color = LightSteelBlue
        t.stroke_line(x, ymin, x, y1)
        t.save_legend_info('50\% $ M_{tot} $')
        x = xs[d.sx_M.where_closest(d.star_Mass * 0.95)]
        t.stroke_color = MediumSlateBlue
        t.stroke_line(x, ymin, x, y1)
        t.save_legend_info('95\%')
        x = xs[d.sx_M.where_closest(d.star_Mass * 0.999)]
        t.stroke_color = MediumBlue
        t.stroke_line(x, ymin, x, y1)
        t.save_legend_info('99.9\%')
    end
    
    def abundances(xs, title, xlabel, xleft, xright, ymin = -5, ymax = 0.1)
        title = get_star_age_title('Abundances --- Age %s') if title == nil
        show_box_labels(title, xlabel, 'log mass frac')
        t.show_plot('left_boundary' => xleft, 'right_boundary' => xright,
            'top_boundary' => ymax, 'bottom_boundary' => ymin) do
            show_mesh(xs)
            t.show_polyline(xs, d.sx_logXH, BrightBlue, '\it Hydrogen', Line_Type_Solid)
            t.show_polyline(xs, d.sx_logXHE, Goldenrod, '\it Helium', Line_Type_Solid)
            t.show_polyline(xs, d.sx_logXC, Coral, '\it Carbon', Line_Type_Solid)
            t.show_polyline(xs, d.sx_logXN, Lilac, '\it Nitrogen', Line_Type_Solid)
            t.show_polyline(xs, d.sx_logXO, FireBrick, '\it Oxygen', Line_Type_Solid)
            t.show_polyline(xs, d.sx_logXNE, RoyalPurple, '\it Neon', Line_Type_Solid)
            show_mass_locs(xs, ymin) unless xs == d.sx_M
        end
    end
    
    def density(xs, title, xlabel, xleft, xright, ymin = -10, ymax = 10)
        title = get_star_age_title('Density --- Age %s') if title == nil
        show_box_labels(title, xlabel, '\textcolor[rgb]{0.645,0.000,0.129}{log $\rho$}')
        t.show_plot('left_boundary' => xleft, 'right_boundary' => xright,
            'top_boundary' => ymax, 'bottom_boundary' => ymin) do
            show_mesh(xs)
            t.show_polyline(xs, d.sx_RHO.safe_log10, BrickRed, nil, Line_Type_Solid)
            show_mass_locs(xs, ymin) unless xs == d.sx_M
        end
    end
    
    def temperature(xs, title, xlabel, xleft, xright, ymin = 2, ymax = 10)
        title = get_star_age_title('Temperature --- Age %s') if title == nil
        show_box_labels(title, xlabel, '\textcolor[rgb]{0.4,0.0,0.6}{log T}')
        t.show_plot('left_boundary' => xleft, 'right_boundary' => xright,
            'top_boundary' => ymax, 'bottom_boundary' => ymin) do
            #show_mesh(xs)
            t.show_polyline(xs, d.sx_T.safe_log10, RoyalPurple, nil, Line_Type_Dash)
            show_mass_locs(xs, ymin) unless xs == d.sx_M
        end
    end
    
    def rho_T_2Ys(xs, title, xlabel, xleft, xright)
        t.subplot {
            t.yaxis_loc = t.ylabel_side = LEFT;
            t.right_edge_type = AXIS_HIDDEN; density(xs, title, xlabel, xleft, xright) }
        t.subplot {
            t.yaxis_loc = t.ylabel_side = RIGHT;
            t.left_edge_type = AXIS_HIDDEN; temperature(xs, title, xlabel, xleft, xright) }
    end
    
    def convection(xs, title, xlabel, xleft, xright, ymin = -0.3, ymax = 1.5)
        title = get_star_age_title('Convection --- Age %s') if title == nil
        show_box_labels(title, xlabel, '$ dlnT/dlnP $')
        t.show_plot('left_boundary' => xleft, 'right_boundary' => xright,
            'top_boundary' => ymax, 'bottom_boundary' => ymin) do
            show_mesh(xs)
            t.show_polyline(xs, d.sx_GRAD_AD, BrightBlue, '$ \nabla_{ad} $', Line_Type_Solid)
            t.show_polyline(xs, d.sx_GRAD_RAD, Goldenrod, '$ \nabla_{rad} $', Line_Type_Solid)
            t.show_polyline(xs, d.sx_GRAD_STAR, Crimson, '$ \nabla_{actual} $', Line_Type_Dot)
            show_mass_locs(xs, ymin)
        end
    end
    
    def grad_star(xs, title, xlabel, xleft, xright, ymin = -0.01, ymax = 0.5)
        title = get_star_age_title('$ \nabla_{actual} $ --- Age %s') if title == nil
        show_box_labels(title, xlabel, '$ \nabla_{actual} $')
        t.show_plot('left_boundary' => xleft, 'right_boundary' => xright,
            'top_boundary' => ymax, 'bottom_boundary' => ymin) do
            show_mesh(xs)
            t.show_polyline(xs, d.sx_GRAD_STAR, BrightBlue)
            show_mass_locs(xs, ymin)
        end
    end
    
    def ionization_ratios(starting_xs, val, divisor, color, text, type)
        xs = @temp_vec1
        ys = @temp_vec2
        frac_limit = 1e-20
        j = 0
        val.each2_with_index(divisor) do |v,d,i|
            if d > frac_limit
                tmp = v/d
                tmp = 0 if tmp < 0
                tmp = 1 if tmp > 1
                xs[j] = starting_xs[i]
                ys[j] = tmp
                j += 1
            end
        end
        xs.resize(j)
        ys.resize(j)
        t.show_polyline(xs, ys, color, text, type)
    end
    
    def ionization(xs, title, xlabel, xleft, xright, ymin = 0.02, ymax = 1.02)
        title = get_star_age_title('Ionization --- Age %s') if title == nil
        show_box_labels(title, xlabel, 'ionization fraction')
        t.show_plot('left_boundary' => xleft, 'right_boundary' => xright,
            'top_boundary' => ymax, 'bottom_boundary' => ymin) do
            show_mesh(xs)
            ionization_ratios(xs, d.sx_XH2, d.sx_XH, Goldenrod, '$ H_{2} $', Line_Type_Solid)
            ionization_ratios(xs, d.sx_XH0, d.sx_XH, BrightBlue, '$ H^{0} $', Line_Type_Solid)
            ionization_ratios(xs, d.sx_XH_plus, d.sx_XH, Coral, '$ H^{+} $', Line_Type_Solid)
            ionization_ratios(xs, d.sx_XHE0, d.sx_XHE, Lilac, '$ He^{0} $', Line_Type_Solid)
            ionization_ratios(xs, d.sx_XHE_plus1, d.sx_XHE, FireBrick, '$ He^{+} $', Line_Type_Solid)
            ionization_ratios(xs, d.sx_XHE_plus2, d.sx_XHE, RoyalPurple, '$ He^{+\!\!+} $', Line_Type_Solid)
            show_mass_locs(xs, ymin)
        end
    end
    
    def ratios_line(xs, vals, divisor, color, text, type)
        ys = @temp_vec2
        ys.replace(vals).div!(divisor).safe_log10!
        t.show_polyline(xs, ys, color, text, type)
    end
    
    def ratios(xs, title, xlabel, xleft, xright, ymin = -5.0, ymax = 0.1)
        title = get_star_age_title('Ratios --- Age %s') if title == nil
        show_box_labels(title, xlabel, 'log ratio')
        t.show_plot('left_boundary' => xleft, 'right_boundary' => xright,
            'top_boundary' => ymax, 'bottom_boundary' => ymin) do
            show_mesh(xs)
            ratios_line(xs, d.sx_T, d.temp_max, BrightBlue, '$ T/T_{max} $', Line_Type_Solid)
            ratios_line(xs, d.sx_OPACITY, d.opacity_max, Goldenrod, '$ \kappa/\kappa_{max} $', Line_Type_Solid)
            ratios_line(xs, d.sx_L, d.luminosity_max, FireBrick, '$ L/L_{max} $', Line_Type_Solid)
            ratios_line(xs, d.sx_RHO, d.density_max, Lilac, '$ \rho/\rho_{max} $', Line_Type_Solid)
            ratios_line(xs, d.sx_R, d.sx_R[d.surface_shell], Coral, '$ R/R_{phot} $', Line_Type_Solid)
            if xs == d.sx_logP
                ratios_line(xs, d.sx_M, d.sx_M[d.surface_shell], RoyalPurple, '$ M/M_{tot} $', Line_Type_Solid)
            else
                ratios_line(xs, d.sx_P, d.sx_P[d.center_shell], RoyalPurple, '$ P/P_{max} $', Line_Type_Solid)
            end
            show_mass_locs(xs, ymin)
        end
    end
    
    def power_line(xs, vals, color, text, type)
        ys = @temp_vec2
        ys.replace(vals).safe_log10!
        t.show_polyline(xs, ys, color, text, type)
    end
    
    def power(xs, title, xlabel, xleft, xright, ymin = -6, ymax = 9)
        title = get_star_age_title('Power --- Age %s') if title == nil
        show_box_labels(title, xlabel, 'Energy input ($log$ $ergs/g/s$)')
        t.show_plot('left_boundary' => xleft, 'right_boundary' => xright,
            'top_boundary' => ymax, 'bottom_boundary' => ymin) do
            show_mesh(xs)
            power_line(xs, d.wimp_E_annihilation, Black, '\it WIMP annihilation', Line_Type_Solid)
            power_line(xs, d.sx_EPS_PP, RoyalPurple, '\it PP', Line_Type_Solid)
            power_line(xs, d.sx_EPS_CNO, Goldenrod, '\it CNO', Line_Type_Solid)
            power_line(xs, d.sx_EPS_3A, Coral, '$ 3\alpha $', Line_Type_Solid)
            power_line(xs, d.sx_EPS_AC, SkyBlue, '$ C+\alpha $', Line_Type_Solid)
            power_line(xs, d.sx_EPS_AN, FireBrick, '$ N+\alpha $', Line_Type_Solid)
            power_line(xs, d.sx_EPS_AO, Lilac, '$ O+\alpha $', Line_Type_Solid)
            power_line(xs, d.sx_EPS_CCA, BrightBlue, '$ C+C $', Line_Type_Dot)
            power_line(xs, d.sx_EPS_NEU, BrightBlue, '\it Neutrinos', Line_Type_Solid)
        end
    end
    
    def show_mesh(xs)
        return unless @show_mesh_flag == true
        ymin = t.bounds_ymin; ymax = t.bounds_ymax
        t.line_color = LightGray; lw = t.line_width; t.line_width = 0.6
        xs.each_with_index { |x, i| t.stroke_line(x, ymin, x, ymax) if i.mod(3) == 1 }
        t.line_width = lw
    end
    
    def show_properites(cmd, xs, title, xlabel, xleft, xright, with_legend)
        t.landscape unless t.in_subplot
        show_model_number
        t.top_edge_type = AXIS_LINE_ONLY
        if (with_legend)
            t.rescale(@legend_scale)
            t.subplot('right_margin' => @plot_with_legend_right_margin) do
                background; cmd.call(xs, title, xlabel, xleft, xright)
            end
            t.set_subframe('left_margin' => @legend_left_margin, 'top_margin' => 0.05)
            t.show_legend
        else
            cmd.call(xs, title, xlabel, xleft, xright)
        end
    end
    
    def show_by_mass(cmd, title, with_legend)
        read_data
        xs = d.sx_M
        xmax = xs.max
        xmax = d.initial_Mass
        show_properites(cmd, xs, title, @mass_xlabel, xs.min, xmax, with_legend)
    end

    def show_by_log_mass(cmd, title, with_legend)
        read_data
        xs = d.sx_M.safe_log10
        xmax = xs.max
        show_properites(cmd, xs, title, @log_mass_xlabel, xs.min, xmax, with_legend)
    end

    
    def show_by_logP(cmd, title, with_legend)
        read_data
        xs = d.sx_logP
        show_properites(cmd, xs, title, @logP_xlabel, xs.max, xs.min, with_legend)
    end
    
    def show_by_both(cmd, title_string, ylabel)
        t.landscape
        read_data
        show_model_number
        t.rescale(@legend_scale)
        t.context do
            t.set_subframe('right_margin' => @plot_with_legend_right_margin)
            show_box_labels(get_star_age_title(title_string), nil, ylabel)
            t.subplot('right_margin' => @side_by_side_margin) do
                t.yaxis_loc = LEFT; t.right_edge_type = AXIS_LINE_ONLY
                background; show_by_mass(cmd, nil, false)
            end
            t.reset_legend_info
            t.subplot('left_margin' => @side_by_side_margin) do
                t.yaxis_loc = RIGHT; t.yaxis_type = AXIS_WITH_TICKS_ONLY; t.left_edge_type = AXIS_LINE_ONLY
                background; show_by_logP(cmd, nil, false)
            end
        end
        t.set_subframe('left_margin' => @legend_left_margin, 'top_margin' => 0.05)
        t.show_legend
    end
    
    def show_abundances(xs, title, xlabel, xleft, xright, with_legend = true)
        show_properites(@abundances_plot, xs, title, xlabel, xleft, xright, with_legend)
    end
    
    def abundances_by_mass
        show_by_mass(@abundances_plot, nil, true)
    end

    def abundances_by_logP
        show_by_logP(@abundances_plot, nil, true)
    end
    
    def abundances_by_both
        show_by_both(@abundances_plot, 'Abundances --- Age %s', 'log mass fraction')
    end

    def convection_by_mass
        show_by_mass(@convection_plot, nil, true)
    end
    
    def grad_star_by_mass
        show_by_mass(@grad_star_plot, nil, true)
    end
    
    def grad_star_by_logP
        show_by_logP(@grad_star_plot, nil, true)
    end
    
    def check_convection
        diff = d.sx_GRAD_AD - d.sx_GRAD_RAD
        i = diff.where_first_gt(0)
        if i >= 0 and d.sx_logP[i] < 8
            j = i
            while diff[j] > 0
                j += 1
            end
            j -= 1 if diff[j] < 0
        end
    end

    def convection_by_logP
        show_by_logP(@convection_plot, nil, true)
        check_convection
    end

    def convection_by_both
        show_by_both(@convection_plot, 'Convection --- Age %s', '$ dlnT/dlnP $')
    end

    def ionization_by_logP
        show_by_logP(@ionization_plot, nil, true)
    end

    def ratios_by_mass
        show_by_mass(@ratios_plot, nil, true)
    end

    def ratios_by_logP
        show_by_logP(@ratios_plot, nil, true)
    end

    def power_by_mass
        show_by_log_mass(@power_plot, nil, true)
    end

    def power_by_logP
        show_by_logP(@power_plot, nil, true)
    end
    
    def trio(mass_flag, title, plot1, plot2, plot3)
        read_data
        show_model_number
        t.rescale(@subplot_scale)
        xlabel = (mass_flag)? @mass_xlabel : @logP_xlabel
        num_plots = 3
        t.legend_scale = 0.6
        t.legend_text_ystart = 0.95
        t.legend_text_dy = 1.2
        row_margin = 0.02
        title = get_star_age_title(title)
        t.subplot(t.row_margins('num_rows' => num_plots, 'row' => 1, 'row_margin' => row_margin)) do
            t.xlabel_visible = false
            t.xaxis_type = AXIS_WITH_TICKS_ONLY
            if mass_flag
                show_by_mass(plot1, title, true)
            else
                show_by_logP(plot1, title, true)
            end
        end
        t.reset_legend_info
        t.top_edge_type = AXIS_LINE_ONLY
        t.subplot(t.row_margins('num_rows' => num_plots, 'row' => 2, 'row_margin' => row_margin)) do 
            t.title_visible = t.xlabel_visible = false
            t.xaxis_type = AXIS_WITH_TICKS_ONLY
            if mass_flag
                show_by_mass(plot2, nil, true)
            else
                show_by_logP(plot2, nil, true)
            end
        end
        t.reset_legend_info
        t.subplot(t.row_margins('num_rows' => num_plots, 'row' => 3, 'row_margin' => row_margin)) do 
            t.title_visible = false
            if mass_flag
                show_by_mass(plot3, nil, true)
            else
                show_by_logP(plot3, nil, true)
            end
        end
    end

    def duo(mass_flag, title, plot1, plot2)
        read_data
        show_model_number
        t.rescale(@subplot_scale)
        xlabel = (mass_flag)? @mass_xlabel : @logP_xlabel
        num_plots = 2
        t.legend_scale = 0.6
        t.legend_text_ystart = 0.95
        t.legend_text_dy = 1.2
        row_margin = 0.02
        title = get_star_age_title(title)
        t.subplot(t.row_margins('num_rows' => num_plots, 'row' => 1, 'row_margin' => row_margin)) do
            t.xlabel_visible = false
            t.xaxis_type = AXIS_WITH_TICKS_ONLY
            if mass_flag
                show_by_mass(plot1, title, true)
            else
                show_by_logP(plot1, title, true)
            end
        end
        t.reset_legend_info
        t.subplot(t.row_margins('num_rows' => num_plots, 'row' => 2, 'row_margin' => row_margin)) do 
            t.title_visible = false
            if mass_flag
                show_by_mass(plot2, nil, true)
            else
                show_by_logP(plot2, nil, true)
            end
        end
    end

    def trio_by_mass
        trio(true, 'Profiles --- Age %s',
            @power_plot, @abundances_plot, @ratios_plot)
    end
    
    def trio2_by_mass
        t.landscape
        trio(true, 'Profiles',
            @power_plot, @abundances_plot, @rho_T_plot)
    end

    def trio2_and_history
        plot_and_history { trio2_by_mass }
    end

    def duo_by_mass
        read_data
        duo(true, 'Rates and Abundances --- Age %s',
            @power_plot, @abundances_plot)
    end

    def trio2_by_logP
        trio(false, 'Profiles ($ \log P $) --- Age %s',
            @power_plot, @abundances_plot, @ratios_plot)
    end

    def trio_by_logP
        trio(false, 'Profiles ($ \log P $) --- Age %s',
            @convection_plot, @ionization_plot, @ratios_plot)
    end

    def trio_by_both
        read_data
        t.landscape unless t.in_subplot
        show_model_number
        column_margin = 0.28
        t.rescale(0.75)
        show_box_labels(get_star_age_title('Summary -- Age %s'))
        t.subplot(t.column_margins('num_columns' => 2, 'column' => 1, 'column_margin' => column_margin)) { trio_by_mass }
        t.subplot(t.column_margins('num_columns' => 2, 'column' => 2, 'column_margin' => column_margin)) { trio_by_logP }
    end
    
    def draw_pair(text, val, x_num, y)
        x_text = x_num - 0.48
        t.show_text('text' => text,
                'x' => x_text, 'y' => y, 'justification' => RIGHT_JUSTIFIED)
        t.show_text('text' => '{\sffamily ' + sprintf('%0.4f',val) + '}',
                'x' => x_num, 'y' => y, 'justification' => RIGHT_JUSTIFIED)
    end
    
    def draw_pair_log(text, val, x_num, y)
        draw_pair(text, val.safe_log10, x_num, y)
    end
    
    def info
        t.rescale(0.7)
        xloc = 1.86
        t.show_text('text' => sprintf('Age %s', d.star_Age_String),
                'side' => TOP, 'position' => xloc, 'scale' => 0.9, 'shift' => 2,
                'justification' => RIGHT_JUSTIFIED, 'color' => Black)
        t.show_text('text' => sprintf('Mass %0.2f ($ \mathrm{M_{\odot}} $)', d.star_Mass),
                'side' => TOP, 'position' => xloc, 'scale' => 0.9, 'shift' => 0.75,
                'justification' => RIGHT_JUSTIFIED, 'color' => Black)
        t.show_text('text' => sprintf('Model %i', d.model_Number),
                'side' => BOTTOM, 'position' => xloc, 'scale' => 0.9, 'shift' => 1.75,
                'justification' => RIGHT_JUSTIFIED, 'color' => Black)
        t.rescale(0.58)
        y = 0.97; dy = -t.default_text_height_dy * 1.2
        draw_pair('$ \log L $', d.log_Luminosity, xloc, y); y += dy;
        draw_pair('$ \log T_{eff} $', d.log_surface_Temp, xloc, y); y += dy;
        draw_pair('$ \log T_{c} $', d.log_center_Temp, xloc, y); y += dy;
        draw_pair('$ \log \rho_{c} $', d.log_center_Density, xloc, y); y += dy;
        draw_pair('$ \Psi_{c} $', d.center_Degeneracy, xloc, y); y += dy;
        draw_pair_log('$ \log scp_{c} $', d.sx_SCP[d.center_shell], xloc, y); y += dy;
        draw_pair('He Core', d.mass_He_Core, xloc, y); y += dy;
        draw_pair('C/O Core', d.mass_C_Core, xloc, y); y += dy;
        y += dy/2;
        draw_pair('$ T_{max} $ Mass', d.sx_M[d.shell_with_max_temp], xloc, y); y += dy;
        draw_pair_log('$ T_{max} $ $\log \rho $', d.sx_RHO[d.shell_with_max_temp], xloc, y); y += dy;
        draw_pair_log('$ \log T_{max} $', d.temp_max, xloc, y); y += dy;
        draw_pair_log('$ \log \kappa_{max} $', d.opacity_max, xloc, y); y += dy;
        draw_pair_log('$ \log L_{max} $', d.luminosity_max, xloc, y); y += dy;
        draw_pair('$ \log \rho_{max} $', d.log_center_Density, xloc, y); y += dy;
        draw_pair_log('$ \log R_{max} $', d.sx_R[d.surface_shell], xloc, y); y += dy;
        draw_pair_log('$ \log P_{max}  $', d.sx_P[d.center_shell], xloc, y); y += dy;
        y += dy/2;
        draw_pair_log('$ \log L_{\nu}$', d.power_Neutrinos, xloc, y); y += dy;
        draw_pair_log('$ \log L_{PP}$', d.power_PP, xloc, y); y += dy;
        draw_pair_log('$ \log L_{CNO}$', d.power_CNO, xloc, y); y += dy;
        draw_pair_log('$ \log L_{3 \alpha}$', d.power_3_alpha, xloc, y); y += dy;
        draw_pair_log('$ \log L_{C+\alpha}$', d.power_C_alpha, xloc, y); y += dy;
        draw_pair_log('$ \log L_{N+\alpha}$', d.power_N_alpha, xloc, y); y += dy;
        draw_pair_log('$ \log L_{O+\alpha}$', d.power_O_alpha, xloc, y); y += dy;
        draw_pair_log('$ \log L_{\mathrm{other}}$', d.power_Metal_burn + d.power_Ne_alpha, xloc, y); y += dy;
        y += dy/2;
        draw_pair('Center XH', d.center_H, xloc, y); y += dy;
        draw_pair('XHe', d.center_He, xloc, y); y += dy;
        draw_pair('XC', d.center_C, xloc, y); y += dy;
        draw_pair('XN', d.center_N, xloc, y); y += dy;
        draw_pair('XO', d.center_O, xloc, y); y += dy;
        draw_pair('XNe', d.sx_XNE[d.center_shell], xloc, y); y += dy;
        y += dy/2;
        draw_pair_log('$ \log t_{nuclear}$', d.nuclear_Timescale, xloc, y); y += dy;
        draw_pair_log('$ \log t_{thermal}$', d.thermal_Timescale, xloc, y); y += dy;
        seconds_per_year = 3.155692597e7
        draw_pair_log('$ \log t_{dynamic}$', d.dynamic_Timescale/seconds_per_year, xloc, y); y += dy;
        draw_pair_log('$ \log t_{step}$', d.time_Step, xloc, y); y += dy;
        y += dy/2;
    end

    def full_profile
        t.title_shift += 0.5
        t.subplot('right_margin' => 0.17) { trio_by_both }
        t.set_subframe('left_margin' => 0.8)
        info
    end

end
