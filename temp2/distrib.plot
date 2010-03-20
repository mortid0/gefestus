#!/usr/bin/ruby

param = "plotx" # с чем работаем impulse, conc, temp, ...
args = "" # параметры графика, для стрелок - with vectors, графики - with lines, поля pm3d
comp = "0" # номер компонента с которым работаем 0...nc
xrange = "xrange [-0.2E3:0.2E3]" # размер оси х
yrange = "yrange [0:5.5E-12]" # размер оси у
plotsize = "size 1024,600" # размер картинки на выходе

`mkdir #{param}/pic`

File.open("plot.script","w+") do |file|
#    file.write "set size 1,1\n" # устанавливаем отношение площади графика к площади картинки по х и у
    file.write "set terminal png #{plotsize}\n"
    file.write "set title 'Asteroid counterflow'\n"
#    file.write "set border 3\n set format x ''\n set format y ''\n"
#    file.write "set pm3d map\n set cbrange [0:2E7]\n"
    file.write "set #{xrange}\n"
#    file.write "set #{yrange}\n"
#    file.write "set xtics nomirror\nset ytics nomirror\n" # убираем черточки на осях
#    file.write "unset mx2tics \nunset my2tics\n" # убираем черточки на осях
#    file.write "unset x2tics \nunset y2tics\n" # убираем черточки на осях
#    file.write "unset xtics\n unset ytics\n"
#    file.write "unset xlabel \nunset ylabel\n" # убираем подписи по осям
    240.times do|i|
	file.write "set output '#{param}/pic/#{'0'*(4-i.to_s.size)}#{i}.png'\n"
	num = i*10
       file.write "plot '#{param}/0/#{num}.txt' using 1:2 with lines title '<- Flow', '#{param}/1/#{num}.txt' using 1:2 with lines title '-> Flow', '#{param}/2/#{num}.txt' using 1:2 with lines title 'Clusters' \n"
    end
end
`rm #{param}/pic/*.png`
`gnuplot plot.script`
#`convert -delay 20 #{param}/pic/*.png distrib.avi`