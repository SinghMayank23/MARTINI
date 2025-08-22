(TeX-add-style-hook "hydro"
 (lambda ()
    (LaTeX-add-index-entries
     "Hydrodynamics"
     "background"
     "Hydro data structure"
     "Hydrodynamics!2+1D Hydro!Kolb"
     "Hydrodynamics!2+1D Hydro!Eskola"
     "Hydrodynamics!3+1D Hydro!Nonaka"
     "Hydrodynamics!3+1D Hydro!Schenke")
    (LaTeX-add-labels
     "kolbhydro"
     "eskolahydro"
     "nonakahydro"
     "schenkehydro"
     "structure"
     "kolbdata"
     "eskoladata"
     "nonakadata"
     "schenkedata")))

