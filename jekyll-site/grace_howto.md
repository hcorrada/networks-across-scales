---
layout: default
title: GRACE HOWTO
---

###Logging in

To log into [Grace](http://www.grace.umd.edu/) machines you 
should use an [SSH](http://en.wikipedia.org/wiki/Secure_Shell)
client. (SSH client is a program that runs on your system. For most
\*NIX machines just type ssh in the terminal. In windows you can use
[PuTTY](http://en.wikipedia.org/wiki/PuTTY">PuTTY)). You should
be able to log in using your university [directory ID and password](http://www.oit.umd.edu/password/). There are Solaris and Linux servers accessible
to you but we strongly recommend that you **only use Linux machines**
so for example you should type `ssh yourname@linux.grace.umd.edu`
. 

###Class directory

The class directory is
`/afs/glue.umd.edu/class/fall2014/cmsc/423/0101`. Within the
`student` sub-directory, each of you should already have a
directory set up. This is where you should do all your work. The
`public` directory will contain common files that you might
need for projects or assignments. 

###Editing programs

You can either write your programs locally and then transfer them
to the servers (using [`scp`](http://en.wikipedia.org/wiki/Secure_copy)
or other SSH file transfer tool like [WinSCP](http://en.wikipedia.org/wiki/WinSCP)
for Windows) or just simply edit and compile the programs directly on
the servers. As you need to make sure that your code compiles and
runs on the servers anyway *I recommend that you do everything on
the servers*. You can write the programs with the editor of your
choice but I recommend a simple *NIX editor like [`emacs`](http://en.wikipedia.org/wiki/Emacs),
[`vim`](http://en.wikipedia.org/wiki/Vim_(text_editor)),
`nano` etc.

###Python

We require that you [python](index.html#python) for all programming
projects in class.
There are a number of versions available in grace, but we ask that you use version 2.7
for your work. To run it in grace use command `python2.7`. We will use the `numpy`,
`Bio` and `matplotlib` libraries in class. These are all available in grace. Programming
project descriptions will include any information about extra libraries if needed.

We've installed a couple of other useful packages, namely IPython and the ability to use the IPython
notebook (although it is awfully slow in grace). To access these, run:

{% highlight bash %}
# for bash
source /afs/glue.umd.edu/class/fall2014/cmsc/423/0101/public/python2.7/bin/activate

# for csh
source /afs/glue.umd.edu/class/fall2014/cmsc/423/0101/public/python2.7/bin/activate.csh
{% endhighlight %}

After doing this, the `python` command will use this installation of python2.7. To use the system default again use:

{% highlight bash %}
deactivate
{% endhighlight %}

To start [Ipython](http://ipython.org/) with matplotlib and numpy imported. 

{% highlight bash %}
ipython --pylab=auto
{% endhighlight %}

Although not recommended in grace since it can be very slow (but highly recommended otherwise), you can use the IPython HTML notebook:

{% highlight bash %}
ipython notebook --pylab=inline &

# look for url in messages
firefox http://127.0.0.1:8888/ &
{% endhighlight %}

In order for this to work, you have to make sure you use X11
forwarding. E.g., ssh using

{% highlight bash %}
ssh -Y yourname@linux.grace.umd.edu
{% endhighlight %}

[This command explained](http://explainshell.com/explain?cmd=ssh+-Y+username%40linux.grace.umd.edu)

###Other useful software
Although python is required for all programming projects, there is other
software available. GNU C and C++ (`gcc` and `g++`) compilers are
installed on all these machines. 

Other popular software is already available on these servers. Type
`tap` to see a list of available applications. For example, you
can type `tap java6` to use Sun java compiler and runtime
environment (`javac` and `javac`) and `tap R` to
use the [R](http://en.wikipedia.org/wiki/R_(programming_language))
statistical computing environment. Scripting language interpreters for 
Perl and Ruby are also available.

If you have any questions do not hesitate to send the TA an email
with CMSC423 in the subject.

###More information

* [OIT Knowledge Base](https://www.itsc.umd.edu/). Search for `116375`.
* [OIT Class Cluster info for CS faculty and instructors](http://www.cs.umd.edu/~larry/new-OIT-cluster-info.html)  
* [General Grace online help](http://www.glue.umd.edu/afs/glue.umd.edu/system/info/olh/)  
