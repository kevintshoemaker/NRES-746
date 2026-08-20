Test command: Rscript -e "devtools::test()"

# Project context

This project develops course materials for NRES 746 (graduate-level custom
statistical modeling for ecology/environmental science, UNR), currently
undergoing a full ground-up redesign. The course objective: students learn
to specify and fit **custom statistical models** to observational data in
ecology/environmental science -- probability, likelihood, Bayesian and ML
inference, numerical algorithms, R/Stan implementation, hierarchical
models, and models for non-independent observations (GAMs, time series,
spatial regression).

Students have limited math backgrounds (Calc II is usually the ceiling;
one intro stats course). Materials must be genuinely graduate-level and
push students to the edge of their ability, but cannot presuppose more
math than that. Primary textbook: Bolker, *Ecological Models and Data in
R* (free draft PDF). Materials include R Markdown lecture/lab handouts,
R scripts, and analysis code.

The course runs 15 weeks, meeting Monday and Wednesday for 50-minute
lectures and Tuesday for a 2h45m lab. Each week follows the same
lecture -> lab -> lecture cycle: Monday lecture, Tuesday lab, Wednesday
lecture. A single lecture topic (e.g. one `LECTUREn.Rmd`) commonly spans
both the Monday and Wednesday sessions, with the intervening lab as a
break in the middle -- so the Wednesday WOD and lecture plan may pick up
mid-document rather than starting a fresh topic, and can reference/recap
material from Monday's session.

## About me

- I primarily code in R. I have some experience with Python, Java, Delphi,
  and Fortran, but R is home base.
- My background is in ecology and conservation biology. I'm picking up more
  advanced math/stats theory as I go.

## R style preferences

- **Base R over tidyverse by default.** Only reach for tidyverse functions
  when they're genuinely simpler or more efficient than the base R
  equivalent for that specific task -- not as a default habit.
- **Comment sparingly.** Prefer equations in comments over prose
  explanations where an equation communicates it better.
- **No comment stacking.** Avoid multiple consecutive comment lines, except:
  - directly below a subheading, or
  - directly above a function definition (to describe what it does).
- **Subheadings**: use the convention `# Label ---------------------------------`
  -- the label and the dashes on the *same* line, nothing on the next line.
  Use subheadings sparingly, just enough to make a script easy to jump
  around in Rstudio (not one per code block).

## Course material conventions

These apply to lab handouts, lecture/WOD handouts, and other
student-facing .Rmd documents:

- **Dual output**: Always include both `word_document` and `pdf_document`
  in the YAML, with the inactive one commented out, so I can switch easily.
  For WOD and lab documents, default to `pdf_document` (comment out
  `word_document`). For other student-facing documents (e.g. lecture
  handouts), default to `word_document` (comment out `pdf_document`).
- **Portable formatting only.** Avoid LaTeX-only raw commands
  (`\newpage`, `\vspace{}`, `\underline{\hspace{}}`, etc.) since they
  silently vanish or render badly in Word output. Use underscore lines
  (e.g. `____________________`) for fill-in blanks and writing space, and
  rely on natural pagination instead of forced page breaks.
- **Instructor-notes block**: include an HTML comment
  (`<!-- ... -->`) near the top of each handout with pacing/timing notes
  and dataset provenance, clearly marked for deletion before printing for
  students. For labs, use this block to also flag the intended
  productive-struggle points (where light scaffolding is deliberate, not
  an oversight) so I don't accidentally "fix" them later.
- **Pen-and-paper labs**: lab handouts are completed without laptops.
  Code chunks in these .Rmd files generate figures/tables for the printed
  handout (set `echo = FALSE`) -- students are not writing or running R
  during lab, so don't add chunks meant for them to execute. WOD handouts
  follow the same rule when the WOD is a math problem; WODs that are
  explicitly programming/pseudocoding exercises should present pseudocode
  or fill-in-the-blank R skeletons rather than runnable chunks.

## Lecture period structure

Every M/W 50-minute lecture period follows the same three-beat shape:

1. **WOD ("Workout of the Day"), ~5-10 min.** A short, ungraded,
   self-contained problem set on a math or programming topic -- either
   directly load-bearing for that day's lecture content or usefully
   tangential (e.g., a calculus refresher the week before it's needed, a
   base-R idiom that will matter in an upcoming lab). Students work alone
   for a few minutes, then briefly compare with a neighbor. WODs should
   rotate across topic types (math derivation, probability puzzle, R/
   pseudocode, numerical-algorithm trace, notation drill) rather than
   repeating the same flavor two days running -- see `WOD1_Diagnostic.Rmd`
   for the target format: numbered parts with a short label and dashed
   subheading rule, fill-in blanks, and (where used) a self-rating table
   at the end so students can flag what needs more study. Not every WOD
   needs a self-rating table -- reserve it for diagnostic/review WODs, not
   routine ones.
2. **Lecture, ~30-35 min.** Core content delivery. When developing a
   lecture plan (not just a WOD), structure it as a short sequence of
   sub-topics with: the key idea in one sentence, the minimal formal
   development needed (no more than the Calc-II/intro-stats floor
   requires), one worked ecological example, and an explicit callout of
   the single point students are most likely to get stuck on (analogous to
   flagging the origin of a tricky term in a derivation). Anchor new
   concepts in ecological intuition first, then formalize -- don't open
   with bare notation.
3. **Discussion, remaining time.** A short, open-ended prompt or question
   set meant to surface misconceptions and get students talking before
   they leave -- not a review quiz. Good discussion prompts ask students
   to compare interpretations, defend a choice, or connect the day's
   content back to an ongoing dataset/model thread. This is a lighter-
   weight cousin of a lab's "group discussion" phase, scaled to a few
   minutes.

When asked to build a lecture plan or WOD, produce all three pieces
together (WOD handout, lecture outline, discussion prompt) unless told
otherwise, and note approximate timing for each in an instructor-notes
block.

## Lab structure: Diagnose -> Derive -> Defend

Labs are built as **three handouts** (Act 1: Diagnose, Act 2: Derive, Act
3: Defend), each following the same individual -> group -> individual
cycle:

- **~15 min individual work**
- **~15 min group discussion**
- **~15 min individual revise-and-reflect**

Roughly 45 min per Act, three Acts per lab (fits the 2h45m period with
room for transitions/wrap-up). This structure suits weeks building
*toward a single model*. (A separate "Model Court" format -- not
individual/group/reflect but a structured debate between competing
candidate models -- is used instead for weeks with genuine model-vs-model
competition; don't force Diagnose->Derive->Defend onto those weeks.)

Use `Lab1_Act1_Diagnose.Rmd` as the template for formatting all three
Acts:

- **YAML**: `title: "Lab N, Act X: <Name>"`, a one-line `subtitle` posing
  the week's driving question, `author: "NRES 746 -- name: ____________________    group: ______"`,
  `date: ""`, dual word/pdf output as above.
- **Instructor-notes HTML comment** right after setup, covering: per-part
  pacing that sums to ~45 min, dataset provenance/source citation, and any
  deliberate scaffolding-reduction notes (what's *supposed* to feel
  underdetermined and why, so groups have something real to reconcile at
  discussion time).
- **Opening "scenario" section**: 1-2 short paragraphs of ecological
  framing before any data or math appears, ending on the concrete question
  the Act will address. Data is presented as "real measurements," with a
  `kable()` table immediately following.
- **Numbered, lettered sub-parts** (`## Part 1: ...`, `**1a.**`, `**1b.**`,
  ...) each followed by underscore blank-lines sized to the expected
  answer length (short numeric answers get one blank line; short-essay
  answers get 2-4).
- **Figures**: generated via `echo=FALSE` chunks; blank axes/grids for
  hand-sketching (see the blank-histogram chunk) where the point is for
  the student to draw on the page, not to look at a finished plot.
  `set.seed()` before any jittering so the figure is reproducible.
- **Group discussion section**: a short bulleted list of comparison
  prompts tied directly back to the individual parts above (e.g., "did
  everyone classify the pattern the same way in Part X? If not, what is
  each person pointing to?") -- not generic discussion questions.
- **Revise-and-reflect section**: always ends with (a) a content question
  testing whether the group converged on the right idea, (b) a forward-
  looking "what would you still need to know to do X" question that seeds
  the next Act/lab, and (c) a **muddiest point** question. These closing
  reflection questions are non-negotiable across all three Acts -- they're
  how students self-diagnose what to study harder, so don't drop them for
  space.

## Calibrating difficulty and level

- This is a graduate course: don't round content down to what an intro
  stats course would cover, but don't assume more math than Calc II plus
  one intro stats course either. When in doubt, build the minimum formal
  machinery needed for that week's model and lean on ecological intuition
  to carry the rest.
- Aim for productive struggle, not comfort. Individual-work sections
  should leave some questions genuinely unresolved for an individual
  working alone -- that's what the group-discussion phase is for. If every
  question in an Act is answerable confidently solo, it's pitched too low.
- Always build in a moment for students to compare/defend differing
  answers (peer evaluation) rather than just checking a single correct
  answer -- Model Court weeks make this the whole point; Diagnose->Derive->
  Defend weeks build it into the group-discussion phase of every Act.
- Always close individual reflection with a self-diagnostic question
  (muddiest point / what would you need to know next) so students have an
  explicit, recorded chance to notice their own gaps.

## Before considering a task done

- Actually knit any .Rmd file produced (both output formats where
  feasible) and confirm it renders without errors -- don't hand back
  unverified R Markdown.
- If a chunk generates a plot meant to illustrate a specific pattern
  (e.g., a candidate curve overlaid on real data), actually render and
  look at it before finalizing parameter choices -- don't guess numbers
  that "should" look right.
- Flag any computed values that come out ugly (long decimals, awkward
  fractions) in places where students will do the same calculation by
  hand -- I'd rather adjust the example than hand students messy
  arithmetic.
- For labs, confirm all three Acts stay on one running dataset/narrative
  thread and that the closing reflection questions in Act 3 genuinely tie
  back to the opening scenario in Act 1.
