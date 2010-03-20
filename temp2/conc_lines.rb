#!/usr/bin/ruby

total = 470
step = 10
start = 0
`mkdir conc/sum`
`mkdir conc/sum/pic`


plot_file = File.open("plot.script","w+")
plot_file.print "set size 1,1\n"
plot_file.print "set terminal png size 1024,600\n"
plot_file.print "set title 'Sum conc of counterflow'\nset yrange[0:7E-14]\n"

total.times{|cur_num|
#=begin
	f1 = File.open("conc/0/#{step*cur_num + start}.txt","r")
	a1 = f1.readlines
	f2 = File.open("conc/1/#{step*cur_num + start}.txt","r")
	a2 = f2.readlines
	f3 = File.open("conc/2/#{step*cur_num + start}.txt","r")
	a3 = f3.readlines

	sum1 = Array.new(a1[0].split.size)
	sum2 = Array.new(sum1.size)
	sum3 = Array.new(sum1.size)
	sum1.size.times{|i| 
		sum3[i] = sum2[i] = sum1[i] = 0.0
	}

	a1.size.times{|i|
		d1 = a1[i].split
		d2 = a2[i].split
		d3 = a3[i].split
		d1.size.times{|j|
			sum1[j] += d1[j].to_f
			sum2[j] += d2[j].to_f
			sum3[j] += d3[j].to_f
		}
	}
	File.open("conc/sum/#{cur_num + start/step}.txt","w+") do |f_out|
		sum1.size.times{|i|
			f_out.print "#{i}; #{sum1[i]}; #{sum2[i]}; #{sum3[i]}; #{sum1[i]+sum2[i]+sum3[i]}\n"
		}
	end
#=end
	plot_file.print "set output 'conc/sum/pic/#{cur_num+start/step}.png'\n"
	plot_file.print "plot 'conc/sum/#{cur_num+start/step}.txt' u 1:2 w l lw 2 ti '<-', 'conc/sum/#{cur_num+start/step}.txt' u 1:3 w l lw 2 ti '->','conc/sum/#{cur_num+start/step}.txt' u 1:4 w l lw 2 ti 'Clusters'\n"
}
plot_file.close

`rm conc/sum/pic/*.*`
`gnuplot plot.script`