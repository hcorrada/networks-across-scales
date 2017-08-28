args = commandArgs(TRUE)
title <- args[1]
#rmd <- args[2]

dir.create(title)
dir.create(file.path(title, "img"))
file.copy(file.path("img", "epiviz.png"), file.path(title, "img", "epiviz.png"))
file.copy(file.path("img", "logo.png"), file.path(title, "img", "logo.png"))
file.copy("template.rmd", file.path(title, "index.rmd"))
file.copy("custom.css", file.path(title, "custom.css"))
file.copy("custom.html", file.path(title, "custom.html"))
file.symlink(path.expand("../notes/images"), title)

#cmd <- paste("cat", "template.rmd", rmd, ">", file.path(title, "index.rmd"))
#system(cmd)
