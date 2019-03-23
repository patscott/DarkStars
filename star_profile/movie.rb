# movie.rb

load 'load.rb'   # use load instead of require so can edit and reload

FigureMaker.default.save_dir = 'movie_out'
FigureMaker.default.add_model_number = true
FigureMaker.default.auto_refresh_filename = '../run/EZ_status.log'
ProfilePlots.new('../run/EZ_status.log', 'movie_out')
