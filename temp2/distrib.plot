#!/usr/bin/ruby

param = "plotx" # � ��� �������� impulse, conc, temp, ...
args = "" # ��������� �������, ��� ������� - with vectors, ������� - with lines, ���� pm3d
comp = "0" # ����� ���������� � ������� �������� 0...nc
xrange = "xrange [-0.2E3:0.2E3]" # ������ ��� �
yrange = "yrange [0:5.5E-12]" # ������ ��� �
plotsize = "size 1024,600" # ������ �������� �� ������

`mkdir #{param}/pic`

File.open("plot.script","w+") do |file|
#    file.write "set size 1,1\n" # ������������� ��������� ������� ������� � ������� �������� �� � � �
    file.write "set terminal png #{plotsize}\n"
    file.write "set title 'Asteroid counterflow'\n"
#    file.write "set border 3\n set format x ''\n set format y ''\n"
#    file.write "set pm3d map\n set cbrange [0:2E7]\n"
    file.write "set #{xrange}\n"
#    file.write "set #{yrange}\n"
#    file.write "set xtics nomirror\nset ytics nomirror\n" # ������� �������� �� ����
#    file.write "unset mx2tics \nunset my2tics\n" # ������� �������� �� ����
#    file.write "unset x2tics \nunset y2tics\n" # ������� �������� �� ����
#    file.write "unset xtics\n unset ytics\n"
#    file.write "unset xlabel \nunset ylabel\n" # ������� ������� �� ����
    240.times do|i|
	file.write "set output '#{param}/pic/#{'0'*(4-i.to_s.size)}#{i}.png'\n"
	num = i*10
       file.write "plot '#{param}/0/#{num}.txt' using 1:2 with lines title '<- Flow', '#{param}/1/#{num}.txt' using 1:2 with lines title '-> Flow', '#{param}/2/#{num}.txt' using 1:2 with lines title 'Clusters' \n"
    end
end
`rm #{param}/pic/*.png`
`gnuplot plot.script`
#`convert -delay 20 #{param}/pic/*.png distrib.avi`