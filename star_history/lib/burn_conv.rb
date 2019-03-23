# burn_conv.rb

class StarHistory

    include Math
    include FigureConstants
    
    def draw_conv_info(xs, s1, e1, s2, e2,s3, e3)
        null_zone = -20
        t.stroke_color = BrightBlue
        t.stroke_width = 2
        xs.each_with_index do |x,i|
            y = e1[i]; t.stroke_line(x, s1[i], x, y) if y > null_zone
            y = e2[i]; t.stroke_line(x, s2[i], x, y) if y > null_zone
            y = e3[i]; t.stroke_line(x, s3[i], x, y) if y > null_zone
        end
    end
    
    def plot_conv_by_mass(xs, show_eps_num=false)
        if show_eps_num
            ylabel = 'with burn zones $\mathrm{M/M_{\\odot}}$'
        else
            ylabel = 'convection zones $\mathrm{M/M_{\\odot}}$'
        end
        t.show_ylabel(ylabel); t.no_ylabel
        conv_M_s1 = d.conv_M_s1[@track_first .. @track_last]
        conv_M_e1 = d.conv_M_e1[@track_first .. @track_last]
        conv_M_s2 = d.conv_M_s2[@track_first .. @track_last]
        conv_M_e2 = d.conv_M_e2[@track_first .. @track_last]
        conv_M_s3 = d.conv_M_s3[@track_first .. @track_last]
        conv_M_e3 = d.conv_M_e3[@track_first .. @track_last]
        ary = [ conv_M_s1, conv_M_e1, conv_M_s2, conv_M_e2, conv_M_s3, conv_M_e3 ]
        ymax = max_of_many(ary) * 1.05
        t.show_plot('boundaries' => [ @start_age, @end_age, ymax, 0.001 ]) do
            background
            mark_profiles_on_x
            t.show_xlabel(@age_xlabel); t.no_xlabel
            draw_conv_info(xs, conv_M_s1, conv_M_e1, conv_M_s2, conv_M_e2, conv_M_s3, conv_M_e3)
            plot_eps_nuc(xs, ymax, false) if show_eps_num
            draw_mesh_lines(xs)
        end
    end
    
    def draw_mesh_lines(xs)
        return unless @mesh_lines_flag
        lw = t.line_width; t.line_width = 0.6
        d.star_Mesh.each { |mesh| stroke_track(xs,mesh[@track_first .. @track_last], DarkGray) }
        t.line_width = lw
    end
    
    def draw_eps_nuc_info(xs, burn1, burn2, burn3, burn4)
        null_zone = -20
        t.stroke_color = Blue
        t.stroke_width = 2
        xs.each_with_index do |x,i|
            z1 = burn1[i]; z2 = burn2[i]; z4 = burn4[i]
            if z2 == null_zone
                t.stroke_line(x, z1, x, z4) if z1 > null_zone
            else
                z3 = burn3[i]
                t.stroke_line(x, z1, x, z2) if z1 > null_zone and z1 < z2
                t.stroke_line(x, z3, x, z4) if z4 > z3
                t.stroke_color = Red
                t.stroke_line(x, z2, x, z3)
                t.stroke_color = Blue
            end
        end
    end
    
    def plot_eps_nuc(xs, ymax = nil, do_background = true)
        t.show_ylabel('burn zone $\mathrm{M/M_{\\odot}}$'); t.no_ylabel
        epsnuc_M_1 = d.epsnuc_M_1[@track_first .. @track_last]
        epsnuc_M_2 = d.epsnuc_M_2[@track_first .. @track_last]
        epsnuc_M_3 = d.epsnuc_M_3[@track_first .. @track_last]
        epsnuc_M_4 = d.epsnuc_M_4[@track_first .. @track_last]
        epsnuc_M_5 = d.epsnuc_M_5[@track_first .. @track_last]
        epsnuc_M_6 = d.epsnuc_M_6[@track_first .. @track_last]
        epsnuc_M_7 = d.epsnuc_M_7[@track_first .. @track_last]
        epsnuc_M_8 = d.epsnuc_M_8[@track_first .. @track_last]
        epsnuc_M_9 = d.epsnuc_M_9[@track_first .. @track_last]
        epsnuc_M_10 = d.epsnuc_M_10[@track_first .. @track_last]
        epsnuc_M_11 = d.epsnuc_M_11[@track_first .. @track_last]
        epsnuc_M_12 = d.epsnuc_M_12[@track_first .. @track_last]
        ary = [ epsnuc_M_1, epsnuc_M_2, epsnuc_M_3, epsnuc_M_4, epsnuc_M_5, epsnuc_M_6, epsnuc_M_7, epsnuc_M_8,
                    epsnuc_M_9, epsnuc_M_10, epsnuc_M_11, epsnuc_M_12 ]
        ymax = max_of_many(ary) * 1.05 if ymax == nil
        t.show_plot('boundaries' => [ @start_age, @end_age, ymax, 0.001 ]) do
            background if do_background
            mark_profiles_on_x
            draw_eps_nuc_info(xs, epsnuc_M_1, epsnuc_M_2, epsnuc_M_3, epsnuc_M_4)
            draw_eps_nuc_info(xs, epsnuc_M_5, epsnuc_M_6, epsnuc_M_7, epsnuc_M_8)
            draw_eps_nuc_info(xs, epsnuc_M_9, epsnuc_M_10, epsnuc_M_11, epsnuc_M_12)
            draw_mesh_lines(xs)
            t.show_xlabel(@age_xlabel)
        end
    end
    
    def plot_convection_epsnuc
        setup_data
        set_track_for_plots
        t.rescale(0.7)
        t.yaxis_numeric_label_scale = 0.55
        t.ylabel_scale = 0.66
        num_plots = 2; row = 1
        xs = d.star_Age[@track_first .. @track_last]
        t.subplot(t.row_margins('num_rows' => num_plots, 'row' => row)) do 
            t.xaxis_type = AXIS_WITH_TICKS_ONLY
            t.no_xlabel
            plot_conv_by_mass(xs)
        end
        row += 1
        if false
            t.subplot(t.row_margins('num_rows' => num_plots, 'row' => 1)) do # first
                t.xaxis_type = AXIS_WITH_TICKS_ONLY
                 t.top_edge_type = AXIS_HIDDEN
               t.no_xlabel
                plot_eps_nuc(xs)
            end
            row += 1
        end
        t.subplot(t.row_margins('num_rows' => num_plots, 'row' => row)) do 
            t.top_edge_type = AXIS_HIDDEN
            plot_conv_by_mass(xs,true)
        end
        t.show_text('text' => '\textcolor{Blue}{Blue = over 1 erg/g/s}',
            'justification' => LEFT_JUSTIFIED,
            'side' => BOTTOM, 'pos' => 0, 'shift' => 2.9, 'scale' => 0.7)
        t.show_text('text' => '\textcolor{Red}{Red = over 1000 erg/g/s}',
            'justification' => RIGHT_JUSTIFIED,
            'side' => BOTTOM, 'pos' => 1, 'shift' => 2.9, 'scale' => 0.7)
    end
    
end
