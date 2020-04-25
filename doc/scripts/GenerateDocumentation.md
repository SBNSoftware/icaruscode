# Rendering and publishing icaruscode documentation {#ICARUS_Doxygen_README}

| Title:          | Rendering and publishing `icaruscode` documentation  |
| --------------- | ---------------------------------------------------- |
| Type:           | Technical documentation                              |
| Author:         | Gianluca Petrillo (petrillo@slac.stanford.edu)       |
| Created on:     | April 23, 2020                                       |
| Version:        | 1.0                                                  |


This document describes a procedure to render the documentation embedded in the
`icaruscode` source code, in Doxygen format, and to publish the result in the
ICARUS documentation web page.


### Format of this document

This file (`README.md`) is written in a dialect of [Markdown] format that is
compatible with [`pandoc`][pandoc] and with Doxygen. The format is designed to
be readable as plain text: the reader will find all the links to external
resources at the end of the document.
The program `pandoc` allows to render this file in other formats. For example,
to render it into Portable Document Format:
    
    pandoc --toc -o README.pdf README.md
    
and to render it into HyperText Transfer Protocol format:
    
    pandoc --toc -o README.html README.md
    


## Outline of the procedure ####################################################


The documentation written in Doxygen format across `icaruscode` source files
can be rendered by the `doxygen` program in one of several formats.
We render it in HTML format, and then copy the result into an area which is
served by Fermilab web server under http://icarus-exp.fnal.gov.

While the procedure is mostly automated, it is requested that an operator pushes
through all the steps:

1. rendering of the documentation
2. transfer of the rendered content to the web server storage area
3. update of web page linking to the available documentation

Each of these steps is carried out with one script.


### Rendering of the documentation

The rendering is the more complex part of the procedure, and also the most
versatile.
The script to run this part is `updateLocalDocumentation.sh`. This script can
be found in `icaruscode` GIT repository under `doc/scripts` directory.
In order to successfully run it, the additional library scripts `settings.sh`
and `utilities.sh` must also be in the same directory as the main script.
There are two possibilities:

* clone the GIT repository altogether, and find the scripts in
  `icaruscode/doc/scripts`;
* download the three scripts from anywhere.

The script, when running, will clone the repository anyway, so the recommended
approach is to clone it first, e.g.:
    
    git clone http://cdcvs.fnal.gov/projects/icaruscode
    icaruscode/doc/scripts/updateLocalDocumentation.sh
    
Note that the script *must* be executed from outside the `icaruscode` GIT
repository, like in the example above.

Normally this will clone or update the repository, and run Doxygen on it.
Note that **update**: the script `updateLocalDocumentation.sh` **does** mess up
with the GIT repository, and while there is some provision not to lose data,
if that repository has anything valuable it should be backed up or pushed before
running.

The script will use a template Doxygen configuration file (by default it's
in `icaruscode/doc`), replace some parameters dynamically and then produce HTML
output into a directory that is configured in the Doxygen file (usually `html`).
A Doxygen "tag" file is also generated and kept in the output directory.
It also copies the configuration file in that output directory, and sets a small
text file `meta/doxygen` with metadata about the documentation generation.
A log file is created in the `log` directory with the full Doxygen output on
screen.
If run again, it will attempt to update the GIT repository and then regenerate
the documentation, overwriting the previous one. Note that there have been
instances where part of the documentation would not be updated; in that case,
the option `--clean` will remove the previous rendering.

As a special mode, the script can be instructed to check out a specific branch
for the documentation, or even to use whatever is already there (option
`--currentbranch`), in which case the repository is not updated. These options
have been tested much less than the standard running mode, so special care
should be taken when running.

At this point, opening the file `html/index.html` with a local browser will
access the documentation. To publish it in the standard ICARUS web page, the two
following steps are now needed.


### Copying the documentation into the standard ICARUS web page

The copy of the documentation rendered in the previous step into the standard
location for the ICARUS web page is performed via `publish.sh` script.
By this time, the script should be available in the same directory as 
`updateLocalDocumentation.sh`, in `icaruscode/doc/scripts` of the cloned GIT
repository.

The data in the web site is organized as follow:
    
    icaruscode
      v08_50_00
        index.html
        meta/doxygen
        ...
      v08_50_01
        index.html
        meta/doxygen
        ...
      latest -> v08_50_01
    
When the script is run, a new directory is added with as name the version of the
code being documented, and the symbolic link `latest` is updated to point to
that version (which, strictly speaking, _is_ the latest being added).

For this operation to succeed, the user must have permissions to write into the
web content directory; and before that, the node where this script is ran must
have POSIX (e.g. NFS) access to that directory. Enabling these is matter of
Fermilab service desk requests.

If the directory already exists, the script will refuse to proceed. This is
intentional, to prevent a potentially complex script from unintentionally
deleting information on the web site in case of errors or bugs.

At this point, the documentation is published in the web but it is not linked
yet in the index page of the various versions. That is the purpose of the last
step.


### Updating the list of documented versions

The HTML page with the list of available versions is generated anew every time
from a HTML template, by the script `updateVersionList.sh`, sibling of the
other two in the same directory.
The page generated in this way is currently served as
https://icarus-exp.fnal.gov/at_work/software/doc/icaruscode/versionlist.html

This script also requires POSIX access to the web storage area.
The template used for the page is "hand-written" copying some boilerplate
HTML from other pages in the ICARUS web site.
There are probably smarter ways to do that, although I am not sure Fermilab is
aware of them.


### Configuration of the scripts

The script `updateLocalDocumentation.sh` supports a few command line options
(run with `--help` to read about them). The other two scripts currently accept
no option at all. All three scripts share some configuration settings, and these
settings are stored in a file aptly called `settings.sh`.

Some parameters as the path to the web storage area, the name of the experiment,
the branch to check out by default, etc., are encoded in the file.


## TL;DR #######################################################################

    
    [[ -d 'icaruscode' ]] || git clone http://cdcvs.fnal.gov/projects/icaruscode
    icaruscode/doc/scripts/updateLocalDocumentation.sh
    icaruscode/doc/scripts/publish.sh
    icaruscode/doc/scripts/updateVersionList.sh
    


## Change log ##################################################################

Version 1.0 (April 23, 2020, petrillo@slac.stanford.edu)
:   original version


[Markdown]: https://daringfireball.net/projects/markdown

[pandoc]: https://pandoc.org/MANUAL.html

[ICARUSdoc]:
  <https://icarus-exp.fnal.gov/at_work/software/doc/icaruscode/versionlist.html>
