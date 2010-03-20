#!/usr/bin/ruby

240.times{|k|
	f1 = File.open("temp/0/#{k*2}.txt","r")
	f2 = File.open("temp/1/#{k*2}.txt","r")
	f3 = File.open("temp/3/#{k*2}.txt","w")
	a = f1.readlines
	b = f2.readlines
	c = Array.new(a.size)
	a.size.times{|i|
		str1 = a[i].split
		str2 = b[i].split
		f3.print "\n" if 0 == a[i].to_f
		f3.print "#{str1[0]} #{str1[1]} #{str1[2].to_f+str2[2].to_f};\n" if (!a[i].empty?)
#		f3.print " \n" if (''==str1[0])
	}
}
