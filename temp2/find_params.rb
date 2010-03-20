#!/usr/bin/ruby
#using SI units
m_earth = 6.0E24
grav_const = 6.67E-11
r_earth = 6_400_000
dens_earth =  m_earth/ (4 * Math::PI * r_earth**3 / 3.0)

m_ast = 6.0E14
kx = 200
ky = 100
n_of_ast = m_earth / m_ast
r_ast = (3. * m_ast / (4. * Math::PI * dens_earth))**(1.0/3.0)
g_ast = 1_000 # avg velocity of Earth
chi = 2.0 * Math::atan(grav_const * 2.0 * m_ast / (4.0*r_ast * g_ast**2))
astr_ed = 10.0E8
tor_vol = Math::PI **2 * 2 * astr_ed * r_earth**2
num_dens_of_ast = n_of_ast / tor_vol
cell_vol = tor_vol / (kx*ky)
size_xyz = cell_vol**(1.0/3.0)

puts "Earth density = #{dens_earth}"
puts "Asteroid properties:"
puts "Number of ast = #{n_of_ast}"
puts "n_ast = #{num_dens_of_ast}"
puts "m_ast = #{m_ast}"
puts "d_ast = #{2.0*r_ast}"
puts "chi = #{chi}rad = #{chi * 180 / Math::PI}deg"
puts "size = #{size_xyz}"
puts "dtp_coll = #{0.3 / (g_ast * num_dens_of_ast * 4.0 * r_ast*r_ast)}"
puts "dtp_move = #{size_xyz / g_ast}"
puts "lambda = #{1.0/(Math::sqrt(2.0) * num_dens_of_ast * 4.0 * r_ast*r_ast)}"#dlina svob probega
puts "dtp_move_l = #{1.0/(Math::sqrt(2.0) * num_dens_of_ast * 4.0 * r_ast*r_ast * g_ast)}"