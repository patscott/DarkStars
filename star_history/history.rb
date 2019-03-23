# history.rb

#$star_data_path = '../run'

load 'load.rb'

file = File.open('../star_profile/path.txt', 'r')
#file.puts($star_data_path)
$star_data_path = file.gets.chop
$star_prof_path = file.gets.chop
file.close
StarHistory.new($star_data_path, $star_prof_path)

