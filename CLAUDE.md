Test command: Rscript -e "devtools::test()"

# Project context

This project develops course materials for NRES 746 (graduate-level custom
statistical modeling for ecology/environmental science, UNR). Students have
limited math backgrounds (Calc II is the ceiling; one intro stats course).
Materials include R Markdown lab handouts, R scripts, and analysis code.

## About me

- I primarily code in R. I have some experience with Python, Java, Delphi,
  and Fortran, but R is home base.
- My background is in ecology and conservation biology. I'm picking up more
  advanced math/stats theory as I go, so explanations that connect back to
  ecological intuition are more useful to me than pure math-for-math's-sake.

## R style preferences

- **Base R over tidyverse by default.** Only reach for tidyverse functions
  when they're genuinely simpler or more efficient than the base R
  equivalent for that specific task — not as a default habit.
- **Comment sparingly.** Prefer equations in comments over prose
  explanations where an equation communicates it better.
- **No comment stacking.** Avoid multiple consecutive comment lines, except:
  - directly below a subheading, or
  - directly above a function definition (to describe what it does).
- **Subheadings**: use the convention `# Label ---------------------------------`
  — the label and the dashes on the *same* line, nothing on the next line.
  Use subheadings sparingly, just enough to make a script easy to jump
  around in Rstudio (not one per code block).

## Course material conventions

These apply to lab handouts and other student-facing .Rmd documents:

- **Dual output**: YAML should default to `word_document`, with
  `pdf_document` present but commented out, so I can switch easily.
- **Portable formatting only.** Avoid LaTeX-only raw commands
  (`\newpage`, `\vspace{}`, `\underline{\hspace{}}`, etc.) since they
  silently vanish or render badly in Word output. Use underscore lines
  (e.g. `____________________`) for fill-in blanks and writing space, and
  rely on natural pagination instead of forced page breaks.
- **Instructor-notes block**: include an HTML comment
  (`<!-- ... -->`) near the top of each handout with pacing/timing notes
  and dataset provenance, clearly marked for deletion before printing for
  students.
- **Pen-and-paper labs**: lab handouts are completed without laptops.
  Code chunks in these .Rmd files generate figures/tables for the printed
  handout (set `echo = FALSE`) — students are not writing or running R
  during lab, so don't add chunks meant for them to execute.
- **Lab structure**: labs are built as "Acts" (e.g., Diagnose → Derive →
  Defend) on an individual → group discussion → individual revise-and-
  reflect cycle, usually threaded through one running real ecological
  dataset for narrative continuity across acts and weeks.

## Before considering a task done

- Actually knit any .Rmd file produced (both output formats where
  feasible) and confirm it renders without errors — don't hand back
  unverified R Markdown.
- If a chunk generates a plot meant to illustrate a specific pattern
  (e.g., a candidate curve overlaid on real data), actually render and
  look at it before finalizing parameter choices — don't guess numbers
  that "should" look right.
- Flag any computed values that come out ugly (long decimals, awkward
  fractions) in places where students will do the same calculation by
  hand — I'd rather adjust the example than hand students messy
  arithmetic.
