# tracks.rb

class StarHistory

    include Tioga
    include Math
    include FigureConstants
    
    def mark_point(x, y, dx, dy, n)
        mark_spot(x, y)
        t.show_label('text' => sprintf('%i', n), 'x' => x+dx, 'y' => y+dy,
            'color' => Crimson,
            'scale' => 0.7, 'justification' => CENTERED, 'alignment' => ALIGNED_AT_MIDHEIGHT)
    end
    
    def mark_he(xs, ys)
        if @he_burn_start != nil
            dy = 0.6 * t.default_text_height_dy
            if @gap_in_tracks
                g = @helium_gap_start
                t.show_label('text' => 'He Flash', 'x' => xs[g-1], 'y' => ys[g-1]-dy,
                    'color' => Crimson,
                    'scale' => 0.6, 'justification' => CENTERED, 'alignment' => ALIGNED_AT_MIDHEIGHT)
                t.show_label('text' => 'After He Flash', 'x' => xs[g], 'y' => ys[g]-dy,
                    'color' => Crimson,
                    'scale' => 0.6, 'justification' => CENTERED, 'alignment' => ALIGNED_AT_MIDHEIGHT)
            else
                g = @he_burn_start
                t.show_label('text' => 'He Burning', 'x' => xs[g], 'y' => ys[g]-dy,
                    'color' => Crimson,
                    'scale' => 0.6, 'justification' => CENTERED, 'alignment' => ALIGNED_AT_MIDHEIGHT)
            end
        end
    end
    
    def mark_points(xs, ys)
        dy = 0.5 * t.default_text_height_dy
        dx = t.default_text_height_dx * (t.xaxis_reversed ? 0.3 : -0.1)
        mark_point(xs[0], ys[0], dx, dy, 0)
        d.profile_Ages.each_with_index do |age, i|
            mod_num = d.star_Age.where_closest(age)
            mark_point(xs[mod_num], ys[mod_num], dx, dy, i+1)
        end
        mark_he(xs, ys)
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
    
    def track_info
        t.rescale(0.8)
        t.set_bounds('bounds_left' => 0, 'bounds_right' => 1, 'bounds_top' => 1, 'bounds_bottom' => 0)
        y = 0.15; dy = 1.7 * t.default_text_height_dy
        x_num = 0.1; x_age = 1.2; x_mass = 1.9
        t.show_text('text' => '\small ' + @age_xlabel, 'x' => x_age, 'y' => y,
            'justification' => RIGHT_JUSTIFIED)
        t.show_text('text' => '\small Mass', 'x' => x_mass, 'y' => y,
            'justification' => RIGHT_JUSTIFIED)
        y -= dy*0.8
        t.show_text('text' => '0', 'x' => x_num, 'y' => y,
            'justification' => CENTERED)
        t.show_text('text' => sprintf(@age_format, 0.0), 'x' => x_age, 'y' => y,
            'justification' => RIGHT_JUSTIFIED)
        t.show_text('text' => sprintf(@mass_format, d.star_Mass[0]), 'x' => x_mass, 'y' => y,
            'justification' => RIGHT_JUSTIFIED)
        d.profile_Ages.each2_with_index(d.profile_Mass) do |age, mass, i|
            y -= dy
            mod_num = d.star_Age.where_closest(age)
            t.show_text('text' => sprintf('%i', i+1), 'x' => x_num, 'y' => y,
                'justification' => CENTERED)
            t.show_text('text' => sprintf(@age_format, d.star_Age[mod_num]), 'x' => x_age, 'y' => y,
                'justification' => RIGHT_JUSTIFIED)
            t.show_text('text' => sprintf(@mass_format, d.star_Mass[mod_num]), 'x' => x_mass, 'y' => y,
                'justification' => RIGHT_JUSTIFIED)
        end
    end
    
    def dots_in_gap(xs, ys)
        return if !@gap_in_tracks
        t.line_type = Line_Type_Dot
        g = @helium_gap_start
        t.stroke_line(xs[g-1], ys[g-1], xs[g], ys[g])
    end
    
    def plot_H_R(right_margin=0, left_margin=0)
        setup_data
        t.landscape unless t.in_subplot
        t.title_shift += 0.5
        margin = 0.1
        t.subplot('right_margin' => right_margin) do
            background
            xs = d.log_surface_Temp
            tmax = xs.max; tmin = xs.min; xmargin = margin*(tmax-tmin)
            left = tmax + xmargin; right = tmin - xmargin 
            ys = d.log_Luminosity
            lmax = ys.max; lmin = ys.min; ymargin = margin*(lmax-lmin)
            top = lmax + ymargin; bottom = lmin - ymargin
            t.show_plot('left_boundary' => left, 'right_boundary' => right,
                'top_boundary' => top, 'bottom_boundary' => bottom) do
                t.show_xlabel('Surface Temperature (log K)')
                t.show_ylabel('Luminosity (log $L/\mathrm{L_\odot})$')
                add_ZAMS(@zams_log_Surface_Temp, @zams_log_Luminosity)
                stroke_track(xs, ys, Blue, nil, Line_Type_Solid, 0)
                dots_in_gap(xs, ys)
                mark_points(xs, ys)
            end
        end
    end
    
    def plot_T_RHO(right_margin=0, left_margin=0)
        setup_data
        t.landscape unless t.in_subplot
        t.title_shift += 0.5
        margin = 0.1
        t.subplot('right_margin' => right_margin) do
            background
            xs = d.log_center_Density
            dmax = xs.max; dmin = xs.min; xmargin = margin*(dmax-dmin)
            left = dmin - xmargin; right = dmax + xmargin
            ys = d.log_center_Temp
            tmax = ys.max; tmin = ys.min; ymargin = margin*(tmax-tmin)
            top = tmax + ymargin; bottom = tmin - ymargin
            t.show_plot('left_boundary' => left, 'right_boundary' => right,
                'top_boundary' => top, 'bottom_boundary' => bottom) do
                t.show_xlabel('Centre Density (log g\,cm$^{-3}$)')
                t.show_ylabel('Centre Temperature (log K)')
                add_ZAMS(@zams_log_Center_Density, @zams_log_Center_Temp)
                add_PSI(@psi0, 0)
                add_PSI(@psi5, 5)
                add_PSI(@psi100, 100)
                stroke_track(xs, ys, Blue, nil, Line_Type_Solid, 0)
                dots_in_gap(xs, ys)
                mark_points(xs, ys)
            end
        end
    end
    
    def plot_H_R_T_RHO
        setup_data
        t.set_portrait
        t.rescale(0.55)
        right_margin = 0.2
        left_margin = 0.85
        t.subplot('bottom_margin' => 0.57) do
            plot_H_R(right_margin, left_margin)
            t.subfigure('left_margin' => left_margin) { track_info }
        end
        t.subplot('top_margin' => 0.57) do
            plot_T_RHO(right_margin, left_margin)
        end
    end
    
end
