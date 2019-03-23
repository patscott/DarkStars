class PatPlots

    include Math
    include Tioga
    include FigureConstants
       
    attr_accessor :mass, :wimp_Density, :log_L, :log_T, :r_chi

    def initialize
        patPlots_init
        return self
    end
    
    def t
        @figure_maker
    end

    def colours
        @colours
    end

    def background
        t.fill_color = OldLace
        t.fill_frame
    end

    def patPlots_init
        @figure_maker = FigureMaker.default
        t.def_eval_function { |str| eval(str) }

        t.save_dir = 'PatPlots_out'
        t.def_figure("HR") {plot_H_R}
        t.def_figure("rchi") {plot_r_chi}
        
        @margin = 0.1       
        t.def_enter_page_function { enter_page }
        read_ZAMS('../star_history/star_data')
        
        @colours = {}
        @colours[0.6] = BrightBlue
        @colours[0.8] = DarkOrange
        @colours[1.0] = Burgundy
        @colours[2.0] = Periwinkle
        @colours[4.0] = DarkChocolate
#        colours[10.0] = SeaGreen
    end

    def enter_page
        t.page_setup(11*72/2,8.5*72/2)
        t.set_frame_sides(0.15,0.85,0.85,0.15) # left, right, top, bottom in page coords        
    end

    def plot_H_R

        @HR_log = [ # columns in WIMPZAMS.dat
            @mass = Dvector.new,
            @wimp_Density = Dvector.new,
            @log_L = Dvector.new,
            @log_T = Dvector.new ]
        Dvector.read('WIMPZAMS.dat', @HR_log)     

        t.legend_text_dy = 1.4
        t.ylabel_scale = 0.8
        t.xlabel_scale = 0.8
        t.yaxis_numeric_label_scale = 0.6
        t.yaxis_numeric_label_scale = 0.6
        t.show_ylabel('Luminosity ($\log$ $L$/L$_\odot$)')
        t.show_xlabel('Surface Temperature ($\log$ K)')
        
        tmax = log_T.max; tmin = log_T.min; xmargin = @margin*(tmax-tmin)
        left = tmax + xmargin; right = tmin - xmargin 
        lmax = log_L.max; lmin = log_L.min; ymargin = @margin*(lmax-lmin)
        top = lmax + ymargin; bottom = lmin - ymargin
        t.show_plot_with_legend('legend_scale' => 1.3,'legend_left_margin' => 0.1,'plot_right_margin' => 0,
            'legend_top_margin' => 0.5) do
            t.show_plot('left_boundary' => left, 'right_boundary' => right,
                'top_boundary' => top, 'bottom_boundary' => bottom) do
                background
                mass.each_index do |i| mark_spot(log_T[i], log_L[i], colours[mass[i]]) end
                add_info_line(@zams_log_Surface_Temp, @zams_log_Luminosity)
                colours.sort.reverse.each do |m, col| t.save_legend_info('text' => m.to_s+' M$_\odot$', 'marker' => Bullet, 
                   'marker_scale' => 0.4, 'marker_color'  => col, 'line_width' => -1) end
                t.save_legend_info('text' => 'ZAMS')
            end
        end       
    end

    def plot_r_chi

        @rchi_log = [ # columns in ZAMSrchi.dat
            @mass = Dvector.new,
            @wimp_Density = Dvector.new,
            @r_chi = Dvector.new]
        Dvector.read('ZAMSrchi.dat', @rchi_log)

        r_Sun = 0.69598e11
        log_wimp_Density = wimp_Density.safe_log10
        log_r_chi = (r_chi/r_Sun).safe_log10
        
        t.legend_text_dy = 1.4
        t.ylabel_scale = 0.8
        t.xlabel_scale = 0.8
        t.yaxis_numeric_label_scale = 0.6
        t.yaxis_numeric_label_scale = 0.6
        t.show_ylabel('WIMP characteristic radius ($\log$ $r_\chi$/R$_\odot$)')
        t.show_xlabel('WIMP halo density ($\log$ GeV/cm$^3$)')
        
        rhomax = log_wimp_Density.max; rhomin = log_wimp_Density.min; xmargin = @margin*(rhomax-rhomin)
        left = rhomin - xmargin; right = rhomax + xmargin 
        rmax = log_r_chi.max; rmin = log_r_chi.min; ymargin = @margin*(rmax-rmin)
        top = rmax + ymargin; bottom = rmin - ymargin
        t.show_plot_with_legend('legend_scale' => 1.3,'legend_left_margin' => 0.1,'plot_right_margin' => 0,
            'legend_top_margin' => 0) do
            t.show_plot('left_boundary' => left, 'right_boundary' => right,
                'top_boundary' => top, 'bottom_boundary' => bottom) do
                background
                #mass.each_index do |i| mark_spot(log_wimp_Density[i], log_r_chi[i], colours[mass[i]]) end
                log_wimp_Density[31..60].each_index do |i| mark_spot(log_wimp_Density[i+30], log_r_chi[i+30], colours[2.0]) end
                colours.sort.each do |m, col| 
                    sorted_data = [Array.new, Array.new, Array.new]
                    sorted_data[0] = (0..mass.length-1).find_all{|i| mass[i] == m}
                    sorted_data[0].each_index do |i| 
                        sorted_data[1][i] = log_wimp_Density[sorted_data[0][i]]
                        sorted_data[2][i] = log_r_chi[sorted_data[0][i]]
                    end
                    #t.show_polyline(sorted_data[1], sorted_data[2], col, m.to_s+' M$_\odot$')
                end 
            end
        end       

    end

    def add_info_line(xs, ys)
        t.line_type = Line_Type_Dash
        t.line_color = SlateGray
        t.line_width = 1
        t.append_points_to_path(xs, ys)
        t.stroke
    end

    def mark_spot(x, y, color = Red)
        t.show_marker('x' => x, 'y' => y, 'marker' => Bullet, 'scale' => 0.5, 'color' => color);
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

end

PatPlots.new
