#!/usr/bin/ruby

h = Hash.new
h["file_name"] = "data.txt"
s = ""
File.open("data.txt") do |file|
while line = file.gets
s += line
end
end
#puts "Contents\n"
#puts s

geom = Regexp.new("krp *\t*= *\t*([0-9]+)") 
h["krp"] = geom.match(s)[1].to_i

geom = Regexp.new("kx *\t*= *\t*([0-9]+)") 
h["kx"] = geom.match(s)[1].to_i
geom = Regexp.new("ky *\t*= *\t*([0-9]+)") 
h["ky"] = geom.match(s)[1].to_i
geom = Regexp.new("nc *\t*= *\t*([0-9]+)") 
h["nc"] = geom.match(s)[1].to_i

geom = Regexp.new("nm *\t*= *\t*([0-9]+)") 
h["nm"] = geom.match(s)[1].to_i

geom = Regexp.new("SIZEX *\t*= *\t*([0-9]+\.*[0-9]*)") 
h["SIZEX"] = geom.match(s)[1].to_i
geom = Regexp.new("SIZEY *\t*= *\t*([0-9]+\.*[0-9]*)") 
h["SIZEY"] = geom.match(s)[1].to_i
geom = Regexp.new("SIZEZ *\t*= *\t*([0-9]+\.*[0-9]*)") 
h["SIZEZ"] = geom.match(s)[1].to_i

geom = Regexp.new("dtp *\t*= *\t*([0-9]+\.*[0-9]*)") 
h["dtp"] = geom.match(s)[1].to_i
puts h.inspect

cent_x = h["kx"]*0.5*h["SIZEX"]
cent_y = h["ky"]*0.5*h["SIZEZ"]
cell_vol = h["SIZEX"]*h["SIZEY"]*h["SIZEZ"]
r_0 = h["kx"]*0.3*h["SIZEX"]
single_mass = 1.0							# mass of single particle
diam = 1.0									# diam of single particle
v_z = 0.0										# all particles fly in single plane
w_0 = 0.001										# anglular velocity, all particles are rotationg as a solid body

puts "Origin is (#{cent_x}; #{cent_y}) with R = #{r_0}\n"
File::open("data.in", "w+") do |file|
h["krp"].times{|krp|
	h["kx"].times{|kx|
		h["ky"].times{|ky|
			h["nc"].times{|nc|
				h["nm"].times{|nm|
					x = kx*h["SIZEX"] - cent_x
						y = ky*h["SIZEY"] - cent_y
						r = Math::sqrt(x**2 + y**2)
						mass = 0
						mass = 1.0/(1.0+r/r_0)**2.5 if (r<r_0)
						weight = mass / cell_vol
						phi = Math::atan2(y, x)
						v_x = w_0 * r_0 * Math::cos(phi + Math::PI * 0.25)
						v_y = w_0 * r_0 * Math::sin(phi + Math::PI * 0.25)
						file.puts "#{krp} #{kx} #{ky} #{nc} #{nm} #{weight} #{single_mass} #{diam} #{v_x} #{v_y} #{v_z}"
				}
			}
		}
	}
}
end
