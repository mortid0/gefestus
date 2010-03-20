#!/usr/bin/ruby

param = "temp" # � ��� �������� impulse, conc, temp, ...
args = "" # ��������� �������, ��� ������� - with vectors, ������� - with lines, ���� pm3d
comp = "0" # ����� ���������� � ������� �������� 0...nc
xrange = "xrange [0:159.0E8]" # ������ ��� �
yrange = "yrange [0:79.0E8]" # ������ ��� �
plotsize = "size 1024,600" # ������ �������� �� ������

`mkdir #{param}/pic`

File.open("plot.script","w+") do |file|
#    file.write "set size 1,1\n" # ������������� ��������� ������� ������� � ������� �������� �� � � �
    file.write "set terminal png #{plotsize}\n"
    file.write "set title 'Asteroid counterflow'\n"
    file.write "set border 3\n set format x ''\n set format y ''\n"
    file.write "set pm3d map\n set cbrange [0:2.1E-6]\n"
#    file.write "set #{xrange}\n set #{yrange}\n set cbrange [0:2E5]\n"
    file.write "set xtics nomirror\nset ytics nomirror\n" # ������� �������� �� ����
    file.write "unset mx2tics \nunset my2tics\n" # ������� �������� �� ����
    file.write "unset x2tics \nunset y2tics\n" # ������� �������� �� ����
    file.write "unset xtics\n unset ytics\n"
    file.write "unset xlabel \nunset ylabel\n" # ������� ������� �� ����
    338.times do|i|
	file.write "set output '#{param}/pic/#{'0'*(4-i.to_s.size)}#{i}.png'\n"
       file.write "splot '#{param}/#{comp}/#{i*10}.txt' matrix with pm3d title 'KIN ENRG' \n"#:(sgn($3)*(2/pi*atan(abs($3)))**0.4):(sgn($4)*(2/pi*atan(abs($4))**0.4)) #{args} notitle\n"
    end
end
`rm #{param}/pic/*.png`
`gnuplot plot.script`
#`convert -delay 15 #{param}/pic/*.jpg temp.mpg`