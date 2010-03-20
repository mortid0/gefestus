#!/usr/bin/ruby

param = "conc" # с чем работаем impulse, conc, temp, ...
args = "" # параметры графика, для стрелок - with vectors, графики - with lines, поля pm3d
comp = "0" # номер компонента с которым работаем 0...nc
xrange = "xrange [0:159.0E8]" # размер оси х
yrange = "yrange [0:59.0E8]" # размер оси у
plotsize = "size 1224,800" # размер картинки на выходе

`mkdir #{param}/pic`

File.open("plot.script","w+") do |file|
    file.write "set size 1,1\n" # устанавливаем отношение площади графика к площади картинки по х и у
    file.write "set terminal jpeg #{plotsize}\n set title 'Conc field of asteroid counterflow'\n"
    file.write "set border\n set format x ''\n set format y ''\n"
    file.write "set pm3d map\nset palette rgbformulae -3,-3,-3\n"
#    file.write "set #{xrange}\n set #{yrange}\n"
#    file.write "set cbrange [0:1.4E-15]\n"
#    file.write "unset colorbox\n"
    file.write "set xtics nomirror\nset ytics nomirror\n" # убираем черточки на осях
    file.write "unset mx2tics \nunset my2tics\n" # убираем черточки на осях
    file.write "unset x2tics \nunset y2tics\n" # убираем черточки на осях
    file.write "unset xtics\n unset ytics\n"
    file.write "unset xlabel \nunset ylabel\n" # убираем подписи по осям
    100.times do|i|
	file.write "set output '#{param}/pic/#{'0'*(4-i.to_s.size)}#{i}.jpg'\n"
       file.write "splot '#{param}/#{comp}/#{i*10}.txt' matrix with pm3d notitle \n"
       #file.write "splot '#{param}/2/#{i*10}.txt' matrix with pm3d notitle, '#{param}/0/#{i*10}.txt' matrix with pm3d notitle \n"
    end
end
`rm #{param}/pic/*.jpg`
`gnuplot plot.script`
#`convert -delay 10 #{param}/pic/*.jpg conc.mpg`