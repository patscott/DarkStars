# profiles.rb

class ProfilePlots

    include Math
    include Tioga
    include FigureConstants
    
# EoS plotting methods
    
    def clip_eos_image
        t.move_to_point(t.bounds_left, t.bounds_bottom)
        t.append_point_to_path(t.bounds_left, 4.9)
        t.append_point_to_path(2.5, t.bounds_top)
        t.append_point_to_path(t.bounds_right, t.bounds_top)
        t.append_point_to_path(t.bounds_right, 6.75)
        #t.append_point_to_path(5.8, 7.3)
        t.append_point_to_path(2, 4)
        t.append_point_to_path(2, t.bounds_bottom)
        t.close_path
        t.clip
    end
    
    def add_EoS_PSI(ary, psi)
        xs = ary[0]; ys = ary[1]
        j = ys.where_closest(t.bounds_ymin + 0.08 * t.bounds_height)
        t.line_type = Line_Type_Dash
        t.line_color = SlateGray
        t.line_width = 1
        t.append_points_to_path(xs, ys)
        t.stroke
        t.show_label('text' => sprintf('$\Psi$=%i', psi), 'x' => 1.78, 'y' => 5.98,
            'color' => Black,
            'scale' => 0.52, 'justification' => CENTERED, 'alignment' => ALIGNED_AT_TOP)
    end
    
    def eos_image(title, zs)
        t.show_title(title);
        t.show_xlabel('log $\\rho$');
        t.show_ylabel('log T');
        t.show_plot('boundaries' => [eos.eos_xmin, eos.eos_xmax, eos.eos_ymax, eos.eos_ymin]) do
            background
            add_eos_legend
            clip_eos_image
            t.show_image(
                'll' => [eos.eos_xmin, eos.eos_ymin],
                'lr' => [eos.eos_xmax, eos.eos_ymin], 
                'ul' => [eos.eos_xmin, eos.eos_ymax], 
                'color_space' => t.mellow_colormap, 'data' => @image_data, 'value_mask' => 255,
                'w' => eos.eos_data_xlen, 'h' => eos.eos_data_ylen)
            burning_lines
        end
    end
    
    def burning_lines
        t.line_width = 1.3
        # hydrogen burning
        xs = Dvector[
            2.52017632601677,
            2.19002974183719,
            1.92721353278038,
            1.85754846149542,
            1.89373296858277,
            1.88188288995298,
            1.58382886039391,
            1.22009536184736,
            0.902225556208372,
            0.643645234206801,
            0.439634875583648,
            0.27324347689909]
        ys = Dvector[
            6.63715538518845,
            6.77092643979871,
            6.86998190301536,
            6.9467525662574,
            7.11615939572185,
            7.29614958554737,
            7.37559553326357,
            7.42994604476722,
            7.4820229700758,
            7.52749960544409,
            7.56462180156689,
            7.59390377229902]
        t.show_polyline(xs, ys, Black, nil, Line_Type_Dash)
        t.show_label('text' => 'Hydrogen', 'x' => xs[-1]+0.4, 'y' => ys[-1]+0.1,
            'color' => Black,
            'scale' => 0.6, 'justification' => RIGHT_JUSTIFIED, 'alignment' => ALIGNED_AT_BASELINE)
        # helium burning
        xs = Dvector[1.9, 3.0, 4.0, 5.0, 6.0, 6.3]
        ys = Dvector[8.25, 8.21, 8.14, 8.03, 7.9, 7.85]
        t.show_polyline(xs, ys, Black, nil, Line_Type_Dash)
        t.show_label('text' => 'Helium', 'x' => xs[0]+0.3, 'y' => ys[0]+0.1,
            'color' => Black,
            'scale' => 0.6, 'justification' => RIGHT_JUSTIFIED, 'alignment' => ALIGNED_AT_BASELINE)
        # carbon burning
        xs = Dvector[3.9, 4.5, 5.0, 6.0, 7.0, 8.0]
        ys = Dvector[9.0, 8.9, 8.85, 8.78, 8.74, 8.72]
        t.show_polyline(xs, ys, Black, nil, Line_Type_Dash)
        t.show_label('text' => 'Carbon', 'x' => xs[0]+0.4, 'y' => ys[0]+0.1,
            'color' => Black,
            'scale' => 0.6, 'justification' => RIGHT_JUSTIFIED, 'alignment' => ALIGNED_AT_BASELINE)
        
        
        return
        # eps_nuc > 10^-8
        xs = Dvector[1.68, 0.14, -2.57]
        ys = Dvector[6.08, 6.33, 6.65]
        t.show_polyline(xs, ys, Black, nil, Line_Type_Dash)
        
        # eps_cno > 10^-8
        xs = Dvector[1.43, 1.12, -0.67, -1.73]
        ys = Dvector[6.74, 6.76, 6.85, 6.92]
        t.show_polyline(xs, ys, Gray, nil, Line_Type_Dash)
        
        
    end
    
    def color_bar(ylabel, levels, colors)
        xmin = 0; xmax = 1; xmid = 0.5
        t.rescale(0.8)
        t.xaxis_type = AXIS_LINE_ONLY
        t.xaxis_loc = BOTTOM
        t.top_edge_type = AXIS_LINE_ONLY
        t.yaxis_loc = t.ylabel_side = RIGHT
        t.yaxis_type = AXIS_WITH_TICKS_AND_NUMERIC_LABELS
        t.left_edge_type = AXIS_WITH_TICKS_ONLY
        t.ylabel_shift += 0.5
        t.yaxis_major_tick_length *= 0.5
        t.yaxis_minor_tick_length *= 0.4
        t.show_ylabel(ylabel); t.no_ylabel
        t.show_plot('boundaries' => [xmin, xmax, @image_zmax, @image_zmin]) do
            t.axial_shading(
                'start_point' => [xmid, @image_zmin], 'end_point' => [xmid, @image_zmax], 
                'colormap' => t.mellow_colormap )
            if levels != nil
                t.line_width = 1.5
                clr_num = 0
                num_colors = colors.size
                levels.each do |level|
                    t.line_color = colors[clr_num]
                    clr_num += 1
                    clr_num = 0 if clr_num >= num_colors
                    t.stroke_line(xmin, level, xmax, level)
                end
            end
        end
    end
    
    def show_section(first_i, last_i, state)
        if first_i == last_i
            if last_i == d.sx_logRHO.size
                first_i = last_i - 1
            else
                last_i += 1
            end
        end
        logRHO = d.sx_logRHO[first_i..last_i]
        logT = d.sx_logT[first_i..last_i]
        if state == 0
            t.line_type = Line_Type_Solid
            t.line_width = 4
            t.show_polyline(logRHO, logT, RoyalPurple)
            t.line_width = 2
            t.show_polyline(logRHO, logT, Teal)
        elsif state == 1
            t.line_type = Line_Type_Solid
            t.line_width = 5
            t.show_polyline(logRHO, logT, RoyalPurple)
            t.line_width = 3
            t.show_polyline(logRHO, logT, Gold)
        else
            t.line_type = Line_Type_Solid
            t.line_width = 6
            t.show_polyline(logRHO, logT, RoyalPurple)
            t.line_width = 4
            t.show_polyline(logRHO, logT, Crimson)
        end
    end
    
    def eos_info_pair(text, val, x_num, y)
        x_text = x_num - 1.5
        t.show_text('text' => text,
                'x' => x_text, 'y' => y, 'justification' => RIGHT_JUSTIFIED)
        t.show_text('text' => '{\sffamily ' + sprintf('%0.2f',val) + '}',
                'x' => x_num, 'y' => y, 'justification' => RIGHT_JUSTIFIED)
    end

    def add_eos_info
        t.context do
            t.rescale(0.6)
            xloc = 8.5; y = 4.2; dy = -t.default_text_height_dy * 1.2
            eos_info_pair('$ \log R_{max} $', d.sx_R[d.surface_shell].safe_log10, xloc, y); y += dy;
            eos_info_pair('He Core', d.mass_He_Core, xloc, y); y += dy;
            eos_info_pair('C/O Core', d.mass_C_Core, xloc, y); y += dy;
        end
    end
    
    def add_eos_legend
        add_eos_info
        x1 = -9; x2 = -8; x3 = -7.4
        dy = -0.07
        y1 = 8.7; y2 = y1 + dy
        t.line_type = Line_Type_Solid
        t.line_width = 5
        t.stroke_color = RoyalPurple
        t.stroke_line(x1, y1, x2, y1)
        t.line_width = 3
        t.stroke_color = Gold
        t.stroke_line(x1, y1, x2, y1)
        t.show_text('text' => 'over 1 erg/g/s', 'at' => [x3, y2],
            'scale' => 0.66, 'justification' => LEFT_JUSTIFIED)
        y1 = 9.0; y2 = y1 + dy
        t.line_type = Line_Type_Solid
        t.line_width = 6
        t.stroke_color = RoyalPurple
        t.stroke_line(x1, y1, x2, y1)
        t.line_width = 4
        t.stroke_color = Crimson
        t.stroke_line(x1, y1, x2, y1)
        t.show_text('text' => 'over 1000 erg/g/s', 'at' => [x3, y2],
            'scale' => 0.66, 'justification' => LEFT_JUSTIFIED)
        dy = -0.09
        y1 = 8.2; y2 = y1 + dy
        xdot = 0.5*(x1+x2)
        draw_mass_spot(xdot, y1, @mass_50_color)
        t.show_text('text' => '50\% $ M_{tot} $', 'at' => [x3, y2],
            'scale' => 0.66, 'justification' => LEFT_JUSTIFIED)
        y1 = 7.9; y2 = y1 + dy
        draw_mass_spot(xdot, y1, @mass_95_color)
        t.show_text('text' => '95\%', 'at' => [x3, y2],
            'scale' => 0.66, 'justification' => LEFT_JUSTIFIED)
        y1 = 7.6; y2 = y1 + dy
        draw_mass_spot(xdot, y1, @mass_999_color)
        t.show_text('text' => '99.9\%', 'at' => [x3, y2],
            'scale' => 0.66, 'justification' => LEFT_JUSTIFIED)
    end
    
    def show_star(which_state)
        eps = d.sx_EPS_NUC # nuclear energy generation rate (erg/g/s)
        start_i = 0
        old_state = 0
        eps.each_index do |i|
            rate = eps[i]
            if rate < 1
                new_state = 0
            elsif rate < 1000
                new_state = 1
            else
                new_state = 2
            end
            if i == 0
                old_state = new_state
                start_i = i
            elsif new_state != old_state
                show_section(start_i, i, old_state) if old_state == which_state
                old_state = new_state
                start_i = i
            end
        end
        show_section(start_i, eps.size-1, old_state) if old_state == which_state
    end
    
    def draw_mass_spot(x, y, color)
        t.show_marker('x' => x, 'y' => y, 'mode' => FILL_AND_STROKE, 'stroke_width' => 0.7,
            'marker' => Bullet, 'scale' => 0.8, 'stroke_color' => Black, 'fill_color' => color);
    end
    
    def find0(xx1,yy1,xx2,yy2)
        # returns x where y is 0 on line connecting the points (xx1,yy1) and (xx2,yy2)
        a = (xx1*yy2)-(xx2*yy1)
        b = yy2-yy1
        if ((a.abs > b.abs*1e30) and ((yy1 >= 0 and yy2 <= 0) or (yy1 <= 0 and y2 > 0)))
            xz = 0.5*(xx1+xx2)
        else
            xz = a/b
        end
        return xz
    end
    
    def show_mass_point(frac, color)
        m = d.star_Mass * frac
        i = d.sx_M.where_closest(m)
        j = i+1
        m0 = d.sx_M[i]; m1 = d.sx_M[j]
        if ((m0-m)*(m-m1) < 0)
            j = i-1; m1 = d.sx_M[j]
        end
        rho = find0(d.sx_logRHO[i], m0-m, d.sx_logRHO[j], m1-m)
        temp = find0(d.sx_logT[i], m0-m, d.sx_logT[j], m1-m)
        draw_mass_spot(rho, temp, color)
    end
    
    def eos_contours(zs, levels, colors)
        t.xaxis_type = t.yaxis_type = AXIS_WITH_TICKS_ONLY
        t.no_title; t.no_xlabel; t.no_ylabel
        t.show_plot('boundaries' => [eos.eos_xmin, eos.eos_xmax, eos.eos_ymax, eos.eos_ymin]) do
            clip_eos_image
            t.stroke_color = Gray
            t.line_width = 1
            dest_xs = Dvector.new; dest_ys = Dvector.new; gaps = Array.new
            dict = { 'dest_xs' => dest_xs, 'dest_ys' => dest_ys, 'gaps' => gaps,
                    'xs' => eos.eos_logRHOs, 'ys' => eos.eos_logTs,
                    'data' => zs }
            clr_num = 0
            num_colors = colors.size
            t.line_width = 1.6
            levels.each do |level|
                dict['level'] = level
                t.make_contour(dict)
                t.line_color = colors[clr_num]
                clr_num += 1
                clr_num = 0 if clr_num >= num_colors
                t.append_points_with_gaps_to_path(dest_xs, dest_ys, gaps, true)
                t.stroke
            end
            add_EoS_PSI(eos.psi5, 5)
            show_star(0)
            show_star(1)
            show_star(2)
            show_mass_point(0.5, @mass_50_color)
            show_mass_point(0.95, @mass_95_color)
            show_mass_point(0.999, @mass_999_color)
        end
    end
    
    def eos_plot(title, zs, z_lower_limit, levels, z_upper_limit = nil)
        #t.landscape
        t.set_aspect_ratio(2)
        read_data
        color_bar_label = title
        title = get_star_age_title(title + ' --- Age %s')
        #show_model_number
        @image_zmin = zs.min_gt(z_lower_limit)
        if z_upper_limit != nil
            @image_zmax = zs.max_lt(z_upper_limit)
        else
            @image_zmax = zs.max
        end
        @image_data = t.create_image_data(zs,
            'min_value' => @image_zmin, 'max_value' => @image_zmax, 'masking' => true)
        t.rescale(0.8)
        colors = [ DarkGray, DimGray ]
        t.subplot('right_margin' => 0.07) do
            eos_image(title, zs)
        end
        t.subplot('left_margin' => 0.95, 'top_margin' => 0.05, 'bottom_margin' => 0.05) do
            color_bar(color_bar_label, levels, colors)
        end
        t.subplot('right_margin' => 0.07) do
            eos_contours(zs, levels, colors)
        end
    end

    def pressure_eos
        levels = [ 7.5, 10.5, 13.5, 16.5, 19.5, 22.5 ]
        eos_plot('log P', eos.pressure_EoS, 0.1, levels)
    end

    def energy_eos
        levels = [13.4, 14.7, 16, 18, 22]
        eos_plot('log U', eos.energy_EoS, 0, levels, 21.5)
    end

    def entropy_eos
        levels = [8.8, 9.0, 9.2, 9.4, 9.6]
        eos_plot('log Entropy', eos.entropy_EoS, 6, levels, 11.6)
    end

    def opacity_eos
        levels = [ 0.4, 2.4, 4.4, 5.4]
        eos_plot('log Opacity', eos.opacity_EoS, -5.5, levels)
    end

    def ionization_eos
        levels = [ -10.5, -7.5, -4.5, -2.5, 2.5 ]
        eos_plot('log (bound $e^-$ / free $e^-$)', eos.ionization_EoS, -20, levels)
    end
    
    def psi_eos
        levels = [ -100, -20, -5, 0, 10 ]
        eos_plot('$\Psi$', eos.psi_EoS, -3, levels, 40)
    end
    
    def plot_and_history(&cmd)
        read_hr_data
        t.landscape
        t.rescale(0.7)
        t.show_title(get_star_age_title('Age %s'))
        show_model_number(1, 1.75)
        t.subplot( 'right_margin' => 0.66, 'top_margin' => 0.115, 'bottom_margin' => 0.115 ) { plot_H_R_T_RHO }
        t.subplot( 'left_margin' => 0.45 ) { cmd.call }
    end

    def entropy_eos_and_history
        plot_and_history { entropy_eos }
    end

    def xxxentropy_eos_and_history
        read_hr_data
        t.landscape
        t.rescale(0.7)
        show_model_number(1, 1.75)
        t.subplot( 'right_margin' => 0.66, 'top_margin' => 0.115, 'bottom_margin' => 0.115 ) { plot_H_R_T_RHO }
        t.subplot( 'left_margin' => 0.45 ) { entropy_eos }
    end

    def read_hr_data
        if !t.root_figure || (@have_hr_data && !t.need_to_reload_data)
            return
        end
        t.need_to_reload_data = true
        read_data
        hr.read_history(@hr_name, t.model_number)
        @have_hr_data = true
#        t.model_number = hr.model_Number
#        read_data
        t.need_to_reload_data = false
    end

    def read_PSIs(path = nil)
        path = '' if path == nil
        path = path + '/' if path.length > 0 && path[-1..-1] != '/'
        Dvector.read(path + 'psi0.data', @psi0 = [Dvector.new, Dvector.new], 2)
        Dvector.read(path + 'psi5.data', @psi5 = [Dvector.new, Dvector.new], 2)
        Dvector.read(path + 'psi100.data', @psi100 = [Dvector.new, Dvector.new], 2)
    end
    
    def read_ZAMS(path = nil)
        path = '' if path == nil
        path = path + '/' if path.length > 0 && path[-1..-1] != '/'
        @zams = [
            @zams_log_Center_Density = Dvector.new,
            @zams_log_Center_Temp = Dvector.new,
            @zams_log_Luminosity = Dvector.new,
            @zams_log_Surface_Temp = Dvector.new,
            @zams_Ms = Dvector.new ]
        Dvector.read(path + 'ZAMS.data', @zams, 2)
    end
    
    def stroke_track(xs, ys, color=Black, type=Line_Type_Solid)
        t.stroke_color = color
        t.line_type = type
        t.line_width = 1.7
        t.append_points_to_path(xs,ys)
        t.stroke
        t.show_marker('marker' => Bullet, 'at' => [xs[-1], ys[-1]], 'scale' => 0.6, 'color' => Crimson)
    end
    
    def add_info_line(xs, ys)
        t.line_type = Line_Type_Dash
        t.line_color = SlateGray
        t.line_width = 1
        t.append_points_to_path(xs, ys)
        t.stroke
    end
    
    def add_ZAMS(xs, ys)
        add_info_line(xs, ys)
    end
    
    def add_PSI(ary, psi)
        xs = ary[0]; ys = ary[1]
        add_info_line(xs, ys)
        j = ys.where_closest(t.bounds_ymin + 0.08 * t.bounds_height)
        t.show_label('text' => sprintf('$\Psi$=%i', psi), 'x' => xs[j], 'y' => ys[j],
            'color' => Crimson,
            'scale' => 0.6, 'justification' => CENTERED, 'alignment' => ALIGNED_AT_MIDHEIGHT)
    end
    
    def plot_H_R
        t.title_shift += 0.5
        margin = 0.1
        background
        xs = hr.log_surface_Temp
        tmax = xs.max; tmin = xs.min; xmargin = margin*(tmax-tmin)
        left = tmax + xmargin; right = tmin - xmargin 
        ys = hr.log_Luminosity
        lmax = ys.max; lmin = ys.min; ymargin = margin*(lmax-lmin)
        top = lmax + ymargin; bottom = lmin - ymargin
        t.show_plot('left_boundary' => left, 'right_boundary' => right,
            'top_boundary' => top, 'bottom_boundary' => bottom) do
            t.show_xlabel('log Surface Temperature')
            t.show_ylabel('log Luminosity $\mathrm{L_\odot}$')
            add_ZAMS(@zams_log_Surface_Temp, @zams_log_Luminosity)
            stroke_track(xs, ys, Blue)
        end
    end
    
    def plot_T_RHO
        t.title_shift += 0.5
        margin = 0.1
        background
        xs = hr.log_center_Density
        dmax = xs.max; dmin = xs.min; xmargin = margin*(dmax-dmin)
        left = dmin - xmargin; right = dmax + xmargin
        ys = hr.log_center_Temp
        tmax = ys.max; tmin = ys.min; ymargin = margin*(tmax-tmin)
        top = tmax + ymargin; bottom = tmin - ymargin
        t.show_plot('left_boundary' => left, 'right_boundary' => right,
            'top_boundary' => top, 'bottom_boundary' => bottom) do
            t.show_xlabel('log Center Density')
            t.show_ylabel('log Center Temperature')
            add_ZAMS(@zams_log_Center_Density, @zams_log_Center_Temp)
            add_PSI(@psi0, 0)
            add_PSI(@psi5, 5)
            add_PSI(@psi100, 100)
            stroke_track(xs, ys, Blue)
        end
    end
        
    def plot_H_R_T_RHO
        read_hr_data
        t.set_portrait
        t.rescale(0.55)
        t.subplot('bottom_margin' => 0.57) { plot_H_R }
        t.subplot('top_margin' => 0.57) { plot_T_RHO }
        xloc = 1.0
        return if t.in_subplot
        t.show_text('text' => sprintf('Age %s', hr.star_Age_String),
                'side' => TOP, 'position' => xloc, 'scale' => 0.9, 'shift' => 2.75,
                'justification' => RIGHT_JUSTIFIED, 'color' => Black)
        t.show_text('text' => sprintf('Mass %0.2f ($ \mathrm{M_{\odot}} $)', hr.star_Mass),
                'side' => TOP, 'position' => xloc, 'scale' => 0.9, 'shift' => 1.1,
                'justification' => RIGHT_JUSTIFIED, 'color' => Black)
        t.show_text('text' => sprintf('Model %i', hr.model_Number),
                'side' => BOTTOM, 'position' => xloc, 'scale' => 0.9, 'shift' => 2.25,
                'justification' => RIGHT_JUSTIFIED, 'color' => Black)
    end

end
