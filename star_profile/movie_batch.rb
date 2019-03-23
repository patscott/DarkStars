# movie_batch.rb

      require 'Tioga/FigureMaker'
      require 'load.rb'

      class ProfilePlots

          def run(movie_save_dir, file_prefix, file_suffix, 
                  plot_number = 0, first_model = 1, last_model = 400)
              t.save_dir = movie_save_dir
              n = first_model
              while n <= last_model
                  fname = file_prefix + n.to_s + file_suffix
                  begin 
                      # check to see if the model file exists
                      f = File.open(fname, 'r')
                      f.close
                      # tell ProfilePlots to reload the data
                      @have_data = false
                      @profile_name = fname
                      puts "\n" + fname
                      t.make_pdf(plot_number)
                      pdf_name = t.figure_pdf(plot_number)
                      name = pdf_name[0..-5] # remove the '.pdf'
                      name = append_sequence_number_to_name(name, n)
                      syscmd = 'mv ' + pdf_name + ' ' + name + '.pdf'
                      system(syscmd)
                  rescue
                      # end up here if the File.open failed
                      # just continue to check the next possible model number
                  end
                  n = n + 1
              end
              puts "\n"
          end
          
          def append_sequence_number_to_name(name, n)
            # always use at least 4 digits.  add leading 0's if necessary.
            name += '_'
            if n < 10
              name += '000'
            elsif n < 100
              name += '00'
            elsif n < 1000
              name += '0'
            end
            name += n.to_s
          end
 
      end
      
      system "rm movie_out/*.pdf"
      file = File.open('path.txt', 'r')
      $star_data_path = file.gets.chop
      $star_prof_path = file.gets.chop
      $star_movie_path = file.gets.chop
      file.close
      ProfilePlots.new($star_data_path+'/EZ_status.log').run('movie_out', $star_movie_path+'/model_', '.log', 1, 1, 10000)
