#!/usr/bin/ruby

make = `make`
puts make
if 0==$?
    `./main > log`
end