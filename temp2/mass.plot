#!/usr/bin/ruby

param = "full_mass" # с чем работаем impulse, conc, temp, ...
args = "" # параметры графика, для стрелок - with vectors, графики - with lines, поля pm3d
comp = "0" # номер компонента с которым работаем 0...nc
xrange = "xrange [0:159.0E8]" # размер оси х
yrange = "yrange [0:79.0E8]" # размер оси у
plotsize = "size 1024,600" # размер картинки на выходе

`mkdir #{param}/pic`

File.open("plot.script","w+") do |file|
    file.write "set size 1,1\n" # устанавливаем отношение площади графика к площади картинки по х и у
    file.write "set terminal jpeg #{plotsize}\n"
    file.write "set border 3\n set format x ''\n set format y ''\n"
    file.write "set palette rgbformulae 4,-3,0\n"
    file.write "set pm3d map\n"
    file.write "set #{xrange}\n set #{yrange}\n set cbrange [6E19:3E20]\n"
    file.write "set xtics nomirror\nset ytics nomirror\n" # убираем черточки на осях
    file.write "unset mx2tics \nunset my2tics\n" # убираем черточки на осях
    file.write "unset x2tics \nunset y2tics\n" # убираем черточки на осях
    file.write "unset xtics\n unset ytics\n"
    file.write "unset xlabel \nunset ylabel\n" # убираем подписи по осям
    400.times do|i|
	file.write "set output '#{param}/pic/#{'0'*(4-i.to_s.size)}#{i}.jpg'\n"
       file.write "splot '#{param}/#{comp}/#{i*10}.txt' using 1:2:3 with pm3d title 'TEMP' \n"#:(sgn($3)*(2/pi*atan(abs($3)))**0.4):(sgn($4)*(2/pi*atan(abs($4))**0.4)) #{args} notitle\n"
    end
end
`rm #{param}/pic/*.jpg`
`gnuplot plot.script`
#`convert -delay 15 #{param}/pic/*.jpg mass.mpg`