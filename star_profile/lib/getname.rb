# getname.rb

module GetModelName

def GetModelName.getname(num)
    return nil if num <= 0
    file = File.open("path.txt")
    path = file.gets.chop
    path2 = file.gets.chop
    file.close
    file = File.open(sprintf("%s/profiles.log", path))
    fullname = nil
    num.times do
        model = -1
        line = file.scanf("%s %s %i")
        return nil if line == nil
        age = line[0]
        mass = line[1]
        model = line[2]
        name = sprintf("model_%i.log", model)
        fullname = path2 + '/' + name
    end
    file.close
    return fullname
end

end
