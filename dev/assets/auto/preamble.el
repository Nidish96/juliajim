;; -*- lexical-binding: t; -*-

(TeX-add-style-hook
 "preamble"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("amsmath" "") ("amssymb" "") ("wasysym" "") ("xcolor" "")))
   (TeX-run-style-hooks
    "amsmath"
    "amssymb"
    "wasysym"
    "xcolor"))
 :latex)

