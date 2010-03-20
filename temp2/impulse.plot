#!/usr/bin/ruby

param = "impulse" # � ��� �������� impulse, conc, temp, ...
args = "with vectors" # ��������� �������, ��� ������� - with vectors, ������� - with lines, ���� pm3d
comp = "2" # ����� ���������� � ������� �������� 0...nc
xrange = "xrange [0:159.0E8]" # ������ ��� �
yrange = "yrange [0:79.0E8]" # ������ ��� �
plotsize = "size 1240,800" # ������ �������� �� ������

`mkdir #{param}/pic`

File.open("plot.script","w+") do |file|
#    file.write "set size 1,1\n" # ������������� ��������� ������� ������� � ������� �������� �� � � �
    file.write "set terminal png #{plotsize}\n"
    file.write "set title 'Asteroid counterflow'\n"
    file.write "set border 15\n set format x ''\n set format y ''\n"
#    file.write "set pm3d map\n set cbrange [0:1.1E6]\n"
#    file.write "set #{xrange}\n set #{yrange}\n set cbrange [0:2E5]\n"
    file.write "set xtics nomirror\nset ytics nomirror\n" # ������� �������� �� ����
    file.write "unset mx2tics \nunset my2tics\n" # ������� �������� �� ����
    file.write "unset x2tics \nunset y2tics\n" # ������� �������� �� ����
    file.write "unset xtics\n unset ytics\n"
    file.write "unset xlabel \nunset ylabel\n" # ������� ������� �� ����
    60.times do|i|
	file.write "set output '#{param}/pic/#{'0'*(4-i.to_s.size)}#{i}.png'\n"
       file.write "plot '#{param}/#{comp}/#{i*10}.txt' u 1:2:(cos($3)*$4*62000000):(sin($3)*$4*62000000) every 2:2 lc 0 #{args} notitle\n"
    end
end
`rm #{param}/pic/*.png`
`gnuplot plot.script`
#`convert -delay 15 #{param}/pic/*.jpg temp.mpg`