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

## Idea bank (Google Doc)

I keep a running Google Doc with WOD ideas and notes about the class,
updated frequently:
https://docs.google.com/document/d/1ZF8rx3dUwppj6ZV4sJ4vJ53XgG-0V9BNKyDwp-7W4D4/edit

When crafting WOD or lab handouts, check this doc for relevant ideas and
incorporate them where they fit, rather than relying only on this file's
static conventions.

## Working with the Bolker PDF

`bolker_book/emdbook.pdf` has front matter that isn't counted in the
book's own page numbers: the PDF's internal page index runs **8 pages
ahead** of the printed page number (e.g. printed page 124 is PDF page
132). When pulling a page range with the `pages` parameter, account for
this offset or you'll land on the wrong content -- check the printed
page number in the corner of the returned image before trusting it.

When pulling a bestiary formula from Bolker for course materials,
double-check it against the wider literature if it looks ecologically
surprising. Bolker's own parameterizations occasionally diverge from
the more common version used elsewhere under the same name (e.g. his
Holling type IV uses $x^2$ in the numerator, while the more widely-cited
Andrews/Monod-Haldane version most other sources call "Holling type IV"
uses a plain $x$ and behaves differently as $x \to \infty$). Bolker's
version isn't wrong -- it's what's actually in the assigned reading --
but it's worth flagging the discrepancy to students rather than
presenting it as the only definition in circulation.

## About me

- I primarily code in R. I have some experience with Python, Java, Delphi,
  and Fortran, but R is home base.
- My background is in ecology and conservation biology. I'm picking up more
  advanced math/stats theory as I go.

## Writing style

- **Avoid em-dashes and other tells of AI-generated writing.** Use a
  comma, a period, or a plain hyphen with spaces (` -- `, as used
  throughout this file) instead of an em-dash. Also avoid other common AI
  writing tics: "it's not just X, it's Y" constructions, rule-of-three
  rhetorical lists, throat-clearing openers ("It's worth noting that...",
  "Furthermore,"), excessive hedging ("arguably", "in many ways"), and
  overused intensifiers like "exactly" and "genuinely" -- cut them or
  replace with a plainer word unless they're doing real work in the
  sentence.
  This applies to all student-facing prose (handouts, discussion prompts,
  lecture text) -- I'm upfront with students that AI helped produce these
  materials, but I'd rather the writing not constantly signal that on its
  own, since it can make students take the exercises less seriously.

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

## Cleaning up a lecture .Rmd

When I ask to "clean up" a `LECTUREn.Rmd` file, this means, specifically:

1. **Name every code chunk** descriptively (e.g. `salmon-ztest-canned`,
   not bare `{r}`), so the chunk outline is useful for jumping around in
   RStudio.
2. **Suppress warnings and messages from the knitted output** (e.g.
   `message = FALSE` in the global `knitr::opts_chunk$set()`, so
   `library()` load messages and similar noise don't show up in the
   HTML).
3. **Modernize the R code**: `<-` for assignment (not `=`), `TRUE`/`FALSE`
   (not `T`/`F`), `seq_len(n)`/`seq_along(x)` instead of `1:n` in for-loop
   headers (the `1:n` idiom silently breaks when `n` is 0), vectorize
   where it doesn't obscure a point the lecture is actively teaching
   (e.g. a for-loop shown specifically to illustrate simulation/resampling
   mechanics should generally stay a for-loop -- that's the pedagogical
   point), and consistent, readable variable/function names (snake_case,
   no dots like `p.val`, no names that shadow base R functions like `df`
   or `confint`). Base R by default; `dplyr`/`tidyr` are fine for data
   wrangling and `ggplot2` for plotting, but this is not a request to
   convert base R code to tidyverse style.

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
  (e.g. `____________________`) for short fill-in blanks, and
  rely on natural pagination instead of forced page breaks.
- **Lines vs. blank space.** Reserve ruled lines (underscore blanks, or
  the `blank_lines()` helper -- see `Lab1_Act1.Rmd`'s setup
  chunk) for **written answers**: prose, interpretation, short-essay
  responses. For questions asking students to manipulate or derive an
  equation (algebra, calculus, "show your work" on a formula), give
  plain blank space instead -- no ruled lines. In PDF-only WOD/lab docs
  this means a `blank_space()` chunk (draws nothing, just reserves
  `fig.height`) rather than `blank_lines()`; lines read as "write words
  here" and can make students feel like they should be composing
  sentences rather than working through math.
- **Pen-and-paper by default; laptops open only once ideas are
  finalized.** For weeks that are primarily math/derivation, lab
  handouts are completed entirely without laptops: code chunks
  generate figures/tables for the printed handout (`echo = FALSE`),
  and students aren't writing or running R during lab. For weeks with
  real coding content (writing functions, running an algorithm,
  fitting a model), use a **closed -> closed -> open -> closed**
  rhythm within each Act instead of banning laptops outright:
  individual work and the group check-in that follows happen on paper,
  with any algorithm finalized before anyone opens a laptop; only then
  does the handout say **"Laptops open now"** (a bolded banner), after
  which AI tools, `?help`, and web search are explicitly fair game,
  since the thinking is done and what's left is translating it into
  working code. The Act's closing group-discussion-and-revision
  section goes back to paper. See `Lab2_Act1.Rmd`, `Lab2_Act2.Rmd`, and
  `Lab2_Act3.Rmd` for the target format.
- **"Pseudocode" means plain English, not R syntax.** A pseudocode
  task (written on paper, before laptops open) means step-by-step
  prose describing the algorithm ("for each candidate value of a: ...
  then store ... then find the smallest"), never an R-syntax skeleton
  with blanks -- that reads as code, not algorithmic thinking, and
  short-circuits the skill being practiced. If the task is instead
  translating an already-known formula into R syntax by hand (no real
  algorithm design involved), call it "write the R function/code," not
  "pseudocode."
- **Optional AI-assisted vs. self-coded split**, for the laptop-open
  portion of a coding-heavy Act: let each group (not individuals within
  a group) choose one of two paths -- (A) write an AI prompt and record
  the prompt (never the AI's code) plus results and a reflection on
  whether it matched their own plan, or (B) write the R themselves, no
  AI, and reflect on what was hardest. This keeps AI use available
  without making it mandatory, and keeps the graded artifact the
  student's own reasoning, not AI-authored code.
- WOD handouts follow the pen-and-paper rule when the WOD is a math
  problem; WODs that are explicitly programming/pseudocoding exercises
  should present pseudocode (plain English, per above) or fill-in-
  the-blank R skeletons rather than runnable chunks.

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
otherwise, and note approximate timing for each.

## Lab structure: three Acts

Labs are built as **three handouts** (Act 1, Act 2, Act 3 -- not a
formalized Diagnose/Derive/Defend framework; a subtitle can still use
words like "derive" or "defend" descriptively if it fits that week's
content), each following the same individual -> group -> individual
cycle:

- **~15 min individual work**
- **~5 min group discussion and revision** (merged, not two separate
  phases -- see below; this is the closing section, distinct from any
  earlier laptop-open "group check-in" a coding-heavy Act might also have)

Roughly 20 min per Act for the simple individual -> group-revision
pattern, three Acts per lab (fits the 2h45m period with room for
transitions/wrap-up and a longer Act 2). Acts with a laptop-open coding
phase run longer than this, sometimes 40+ minutes, since they also
include a group check-in before laptops open plus the coding phase
itself -- budget accordingly rather than compressing any of it to fit.
This structure suits weeks building
*toward a single model*. (A separate "Model Court" format -- not
individual/group/reflect but a structured debate between competing
candidate models -- is used instead for weeks with genuine model-vs-model
competition; don't force this individual/group/revise structure onto
those weeks.)

Use `Lab1_Act1.Rmd` as the template for formatting all three
Acts:

- **YAML**: `title: "Lab N, Act X: <Name>"`, a one-line `subtitle` posing
  the week's driving question, `author: "NRES 746 -- name: ____________________    group: ______"`,
  `date: ""`, dual word/pdf output as above.
- **Opening "scenario" section**: 1-2 short paragraphs of ecological
  framing before any data or math appears, ending on the concrete question
  the Act will address. Data is presented as "real measurements," with a
  `kable()` table immediately following.
- **Numbered, lettered sub-parts** (`## Part 1: ...`, `**1a.**`, `**1b.**`,
  ...) each followed by a blank sized to the expected answer length
  (short numeric answers get one inline blank; short-essay answers get
  2-4 lines) -- and by *kind*, per the lines-vs-blank-space rule above:
  written/interpretive answers get lines, equation manipulation or
  "show your work" gets plain blank space.
- **Figures**: generated via `echo=FALSE` chunks; blank axes/grids for
  hand-sketching (see the blank-histogram chunk) where the point is for
  the student to draw on the page, not to look at a finished plot.
  `set.seed()` before any jittering so the figure is reproducible.
- **Group discussion and revision section** (merged, not separate
  "group discussion" + "revise and reflect" phases): keep this
  deliberately open-ended rather than scaffolded. A one-line prompt to
  compare answers with the group and revise as needed is usually enough
  -- optionally naming which part(s) are worth focusing on (e.g. "compare
  your derivative and critical point"), but not a bulleted question tied
  to every individual part. Close with a generic instruction to record
  changes, new insights, and any "muddiest points" or topics to review,
  followed by one open block of blank lines (`blank_lines(6)` or similar)
  for all of it. Don't pre-write specific reflection questions (a content
  check, a forward-looking question, a muddiest-point question) -- let
  differences between group members' answers be what drives the
  discussion and reflection, rather than a fixed checklist. See
  `Lab1_Act1.Rmd`'s "Group discussion and revision" section for
  the target format.

## Calibrating difficulty and level

- This is a graduate course: don't round content down to what an intro
  stats course would cover, but don't assume more math than Calc II plus
  one intro stats course either. When in doubt, build the minimum formal
  machinery needed for that week's model and lean on ecological intuition
  to carry the rest.
- **Minimal scaffolding, genuine first attempts.** Pitch questions so
  they require students to stretch beyond routine application, not just
  plug numbers into a given formula. Resist the urge to pre-break a
  problem into hand-held sub-steps or hint at the right path before
  students have tried one themselves -- it's fine, even desirable, for an
  individual's first attempt to be wrong. Students learn more by working
  out where their own reasoning broke, with help from peers and the
  instructor, than by being walked to a correct answer from the start.
  When in doubt, cut scaffolding rather than add it.
- Aim for productive struggle, not comfort. Individual-work sections
  should leave some questions genuinely unresolved for an individual
  working alone -- that's what the group-discussion phase is for. If every
  question in an Act is answerable confidently solo, it's pitched too low.
- Always build in a moment for students to compare/defend differing
  answers (peer evaluation) rather than just checking a single correct
  answer -- Model Court weeks make this the whole point; regular
  three-Act weeks build it into the group-discussion-and-revision phase
  of every Act.
- Give students an explicit, recorded chance to self-diagnose gaps
  (muddiest points, what they'd still need to know) -- but as a generic
  invitation within the open group-discussion-and-revision block, not as
  a separate mandated question. See the "Lab structure: three Acts"
  section above.

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
  thread and that the closing reflection questions in Act 3 tie back to
  the opening scenario in Act 1.
- Verify statistical or mathematical claims that go into student
  materials numerically before writing them down, not just visually --
  don't assert two SSR conventions give the same fit, or that a model
  comparison is "fair," without actually running it and checking.
- If a PDF you already sent the user won't re-knit ("I can't write on
  file ... .pdf"), it's almost always a file lock from their viewer
  having it open, not a content error. Render a scratch copy elsewhere
  to confirm the content is fine, ask them to close the viewer, then
  re-knit to the real path -- don't assume the Rmd itself is broken.
