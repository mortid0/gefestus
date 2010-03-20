#!/usr/bin/ruby

param = "conc" # � ��� �������� impulse, conc, temp, ...
args = "" # ��������� �������, ��� ������� - with vectors, ������� - with lines, ���� pm3d
comp = "0" # ����� ���������� � ������� �������� 0...nc
xrange = "xrange [0:159.0E8]" # ������ ��� �
yrange = "yrange [0:59.0E8]" # ������ ��� �
plotsize = "size 1224,800" # ������ �������� �� ������

`mkdir #{param}/pic`

File.open("plot.script","w+") do |file|
    file.write "set size 1,1\n" # ������������� ��������� ������� ������� � ������� �������� �� � � �
    file.write "set terminal jpeg #{plotsize}\n set title 'Conc field of asteroid counterflow'\n"
    file.write "set border\n set format x ''\n set format y ''\n"
    file.write "set pm3d map\nset palette rgbformulae -3,-3,-3\n"
#    file.write "set #{xrange}\n set #{yrange}\n"
#    file.write "set cbrange [0:1.4E-15]\n"
#    file.write "unset colorbox\n"
    file.write "set xtics nomirror\nset ytics nomirror\n" # ������� �������� �� ����
    file.write "unset mx2tics \nunset my2tics\n" # ������� �������� �� ����
    file.write "unset x2tics \nunset y2tics\n" # ������� �������� �� ����
    file.write "unset xtics\n unset ytics\n"
    file.write "unset xlabel \nunset ylabel\n" # ������� ������� �� ����
    100.times do|i|
	file.write "set output '#{param}/pic/#{'0'*(4-i.to_s.size)}#{i}.jpg'\n"
       file.write "splot '#{param}/#{comp}/#{i*10}.txt' matrix with pm3d notitle \n"
       #file.write "splot '#{param}/2/#{i*10}.txt' matrix with pm3d notitle, '#{param}/0/#{i*10}.txt' matrix with pm3d notitle \n"
    end
end
`rm #{param}/pic/*.jpg`
`gnuplot plot.script`
#`convert -delay 10 #{param}/pic/*.jpg conc.mpg`