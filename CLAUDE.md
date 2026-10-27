# Project context

Course materials for NRES 746 (UNR), a graduate course in custom
statistical modeling for ecology and environmental science, in the
middle of a ground-up redesign. Students learn to specify and fit
**custom models** to observational data: probability, likelihood,
Bayesian and ML inference, numerical algorithms, R/Stan, hierarchical
models, and models for non-independent data (GAMs, time series,
spatial).

- **Students:** Calc II is usually the ceiling, plus one intro stats
  course. Materials must be graduate-level and push them, but can't
  presuppose more math than that.
- **Textbook:** Bolker, *Ecological Models and Data in R* (free draft,
  `bolker_book/emdbook.pdf`). Other readings are in `pdf/`.
- **Schedule:** 15 weeks. Monday and Wednesday 50-min lectures, Tuesday
  2h45m lab. A topic often spans Monday and Wednesday, with the lab in
  between. The schedule of record is `CourseSchedule.csv`.
- **Branch:** work happens on `gh-pages`.
- **R:** `Rscript` isn't on the bash PATH; use
  `"/c/Program Files/R/R-4.5.1/bin/Rscript.exe"`. Pandoc is 1.17 (no
  Lua filters or fenced divs).

**Idea bank** (check it when building WODs and labs):
https://docs.google.com/document/d/1ZF8rx3dUwppj6ZV4sJ4vJ53XgG-0V9BNKyDwp-7W4D4/edit

**Bolker PDF:** PDF page = printed page + 8. Check the printed number
on the returned page. When pulling a bestiary formula, check it
against the wider literature if it looks surprising; Bolker's
parameterizations sometimes differ from the common version (his
Holling type IV has $x^2$ in the numerator; the Andrews/Monod-Haldane
version most sources use has $x$). Use his version (it's the reading)
but flag the difference to students.

**About me:** R is home base (some Python, Java, Delphi, Fortran).
Background in ecology and conservation biology; picking up more math
and stats theory as I go.

# Writing style

Applies to all student-facing prose (web pages, handouts, keys,
discussion prompts). Students know AI helped make the materials,
but writing that signals it makes them take exercises less seriously.

- **No em-dashes**, and no `--` in knitted `.Rmd` files (pandoc turns
  it into a dash). Use a comma, period, or parentheses.
- **Avoid AI tics:** "it's not just X, it's Y"; rule-of-three lists;
  throat-clearing ("It's worth noting that", "Furthermore"); hedging
  ("arguably", "in many ways"); intensifiers ("exactly", "precisely",
  "genuinely") unless they do real work.
- **Banned phrases:** "pin down" (say "how uncertain," "rule out");
  callbacks framed as "X is [earlier material] applied to Y" (say what
  X does). Plain references ("Monday's pika data") are fine in
  handouts, not on web pages. Also "reveal" ("WOD reveal").
- **Plain headers:** say what the section is ("WOD: bias-variance
  tradeoff", "Big picture"), not a teaser or a cute label. In
  student-facing pages, label in-class prompts "Question:" or
  "Discussion:", not "Pause" (instructor notes and scripts keep `[PAUSE]`).
- **Vary the phrasing.** Don't reuse a sentence frame or label section
  after section (e.g. a "Sticking point:" label on every section); avoid
  catchy antitheses ("penalizes how *many*; these penalize how *big*").
- **Web pages** have extra rules (public audience, no day references):
  see "Weekly web page".
- **Concise:** if a shorter phrasing says the same thing, use it.
- **No specific forward references** ("this comes up in Lab 7"), since
  schedules shift; say "later in the semester" unless you're editing
  that material now.

# R style

- Base R by default; tidyverse only when it's clearly simpler for the
  task (`dplyr`/`tidyr` for wrangling and `ggplot2` are fine).
- Comment sparingly; equations over prose. No stacked comment lines,
  except right below a subheading or right above a function definition.
- Subheadings: `# Label ---------------------------------` on one line.
  Use them sparingly, for jumping around in RStudio.
- **"Clean up" a `LECTUREn.Rmd`** means: name every chunk
  descriptively; suppress messages/warnings globally; modernize (`<-`,
  `TRUE`/`FALSE`, `seq_len()`/`seq_along()` in loop headers, snake_case,
  no names that shadow base R like `df`); vectorize unless a loop is
  the point being taught. Not a request to convert to tidyverse.

# Student-facing handouts (WODs, labs)

- **Output:** YAML lists both `pdf_document` and `word_document`, one
  commented out. WODs and labs default to PDF; other student documents
  default to Word.
- **Portable formatting:** no raw LaTeX (`\newpage`, `\vspace`, etc.)
  in handouts; underscore blanks for short answers; natural pagination.
- **Lines vs. space:** ruled lines (`blank_lines()`, see
  `labs/Lab1_Act1.Rmd`) for written answers; plain space
  (`blank_space()`) for math and derivations. Be generous with math
  space; a WOD or Act may run to several pages. Never shrink answer
  space to save pages.
- **Fill-in tables:** use `blank_table()` (see `labs/Lab1_Act2.Rmd`,
  `labs/Lab2_Act3.Rmd`), never a markdown table (rows can't be padded
  portably). `""` for student blanks; plotmath only for single symbols
  in headers, formulas go in the prose above. `kable()` is fine for
  display-only tables.
- **Figures:** `echo = FALSE`; blank axes/grids where students sketch;
  `set.seed()` before any randomness.
- **Data files:** CSV in `data/`; give the URL
  `https://kevintshoemaker.github.io/NRES-746/data/<file>.csv` and a
  `read.csv()` chunk (URL on its own line). Paste-in vectors as a
  fallback for small data. Push to `gh-pages` and check the URL loads
  before the handout goes out (see `labs/Lab6_Act2.Rmd`).

## Laptops and AI

- **Pen and paper by default.** Math weeks are laptop-free; figures and
  tables are pre-rendered. Using `rbinom()`, `dpois()`, `qnorm()` etc.
  is not a reason to open laptops (print the figure or table instead).
- **Coding weeks: closed -> closed -> open -> closed** within an Act.
  Individual work and the group check-in happen on paper, with the
  algorithm finalized; then a bold **"Laptops open now"** banner
  (`?help` and web search allowed); the closing discussion is back on
  paper. See `labs/Lab2_Act1.Rmd`. A WOD can use the same banner for a
  hands-on R section (see `wod/WOD_Week5_Lecture2_WritingTheLikelihood.Rmd`,
  Part 5).
- **Pseudocode = plain English steps**, never an R skeleton with
  blanks. If the task is translating a known formula into R, call it
  "write the R code." A programming WOD without laptops uses
  pseudocode or a fill-in R skeleton, not runnable chunks.
- **AI labels** (see `labs/Lab5_Act2.Rmd`). Every coding question gets
  one, in italics after the number (`**2a.** *(Recommended: no AI.)*`):
  - *Recommended: no AI*: the week's core skill (writing a likelihood,
    calling `optim()`); `?help`, web, and classmates are fine.
  - *AI optional/encouraged*: results matter, the code isn't the point
    (plot formatting, loops over a grid).

  After each coding question: `AI use (circle one): no AI / AI helped
  (e.g., debugging) / AI wrote the code`. Explain the labels at the
  start of the laptop section. Don't add "using AI never costs
  points." Close with one reflection question (where AI helped and
  what you had to fix; or, if none, what was hardest). Keys treat a
  label/disclosure mismatch as a conversation, not a penalty.
- **Students submit outputs, not code:** filled-in values,
  interpretation, and named figures/tables. Prefer outputs that show
  the student's own choices (e.g. their hand-guessed curve next to the
  fit). Final projects are the exception.

# Lecture periods

Each 50-min lecture has three beats.

1. **WOD (~5-10 min).** Short, ungraded, self-contained. Alone, then
   compare with a neighbor. Rotate types (derivation, probability
   puzzle, R/pseudocode, algorithm trace, notation drill). Format: see
   `wod/WOD_Week1_Lecture1_Diagnostic.Rmd` (self-rating tables only for
   diagnostic/review WODs).
2. **Lecture (~30-35 min).** Open with 2-3 min of big picture (which
   workflow box, what new question we can answer, link to last lecture
   and this week's lab), and return to it at the end. Each sub-topic:
   key idea, minimal formal development, one ecological example, the
   one likely sticking point. Intuition first, notation second.
3. **Discussion (remaining time).** Open-ended prompts that make
   students compare or defend interpretations. Put the numbers the
   question turns on into the prompt. Where it should land goes in the
   instructor notes, not on the page.

Build all three together (with timing) unless told otherwise. The
lecture and discussion live on the week's web page; see "Weekly
materials".

**Calibration: lectures cover about half of what's planned.** Weeks
7-8 each covered 3-4 of 7-8 planned sections. Plan **3-4 sub-topics**
(counting the big-picture opening) with one example carried through;
about one algorithm plus its example per ~15 min. Put extras in a
marked "Bonus" section. A derivation eats its whole slot, so put
the worked numeric example first. Leave slack, and check whether this
week's lab uses anything lecture hasn't covered yet. Update this
paragraph when the instructor reports how a lecture went.

**WODs run about double their planned time** (5 planned = 10 real).

- About 2/3 short in-class WODs (one or two parts; cut parts rather
  than space), 1/3 take-home. Skipping the WOD is fine when the lecture
  needs the time.
- Anything involved or with real coding is a **take-home**, preferably
  Wednesday (due Monday). Take-homes use the AI labels and no
  neighbor/group comparison. Budget ~5-10 min of the next lecture to
  go over it.

**Later in the course:** occasional paper-discussion Wednesdays (pick a
paper tied to the week's methods; reading guide; questions like "what
are the deterministic function and error distribution?" and "is the
model corroborated or just fit?"; give it a week ahead), and full
end-to-end analysis workflows (messy data to publication figures),
which double as final-project models.

**Naming.** `wod/WOD_Week<N>_Lecture<M>_<Topic>.Rmd`, with header
`**NRES 746, Week <N>, Lecture <M>: <Title>**`. `Lecture1`/`Lecture2`
= Monday/Wednesday slot. If the schedule shifts, rename the file and
its header together.

# Weekly materials

Adopted Oct. 27, 2026 (Week 10). The instructor teaches from a single
web page per week, walking through it in lectures and in the lab's
lecture portion. An ordinary week has:

1. **WOD handout(s) with keys** (see "Lecture periods", "Answer keys").
2. **One web page** for the week (below). Replaces the old PDF NOTES
   as the main lecture resource.
3. **One R script** with the page's code (below).
4. **Lab handout(s) with keys** (see "Labs", "Answer keys").
5. **Short bulleted instructor notes, one per class meeting** (Monday,
   Tuesday's lab lecture, Wednesday; below).

**Mine old material first:** legacy root-level `LECTUREn.Rmd` /
`LABn.Rmd` examples were tested in class. Also check the Bolker
section, the idea bank, and this week's WOD/lab so the page calls back
to them rather than repeating or contradicting them. Worked examples
must use different numbers/functions from the lab's graded questions.

## Weekly web page

Root-level `Week<N>_<Topic>.Rmd` (model:
`Week10_ModelComparison.Rmd`), public on the course site.

- **YAML:** `html_document`, `theme: flatly` (the homepage's theme, not
  the legacy `cerulean`), `toc: yes`, `toc_float: yes`, `toc_depth: 3`.
  Render with `rmarkdown::render_site("Week<N>_<Topic>.Rmd")` so it gets
  the site navbar.
- **Navbar:** add the page to the Lectures menu in `_site.yml`, then
  re-render `index.Rmd` too (each page's navbar is baked in at render
  time; other pages pick up the link only when re-rendered).
- **Structure:** a link to the companion script and the relevant Bolker
  section at the top (no other readings; the instructor assigns those
  on Canvas);
  one H1 per class meeting, titled by topic (no day names), each opening
  with a 2-3 sentence summary; numbered H2 sections inside; "Bonus" for
  extras. The lab's lecture portion gets its own H1; lab Act answers and
  anything the Acts ask students to work out stay off the page.
- **Public audience:** no class-day or week references ("Monday's
  creepers", "Week 7's cutoff"); link the relevant section instead,
  with explicit `{#id}` anchors on referenced headers. Explain course
  jargon on first use (e.g. what a WOD is).
- **Prompts:** `> **Question:**` blockquotes and a `## Discussion`
  section, with no answers on the page (answers go in the instructor
  notes). No timings or instructor-only text; there is no instructor
  version of the page.
- **Code:** demo code in visible chunks; figure helpers and
  figure-only code in `echo = FALSE` chunks. Numbers in the prose come
  from inline R.
- **Widgets** (interactive sliders, where they help): standalone
  HTML/JS files in `widgets/` (models: `lasso_kink.html`,
  `ridge_prior.html`), included with
  `htmltools::includeHTML(file.path("widgets", file))`. Plain
  JavaScript, no Shiny (no server on GitHub Pages). Scope CSS under a
  class and prefix element ids so several widgets can share a page.
- **Stan chunks:** `cache = TRUE`, and call `$draws()` on each fit
  inside the same chunk (cached fits otherwise point at deleted temp
  CSVs). Hide that line with `echo = -n`; knitr's expression count can
  be off by one or two, so check the HTML.
- **Check before committing:** view it through the local server
  (`.claude/launch.json` config `site`, port 8746; a `file://` URL
  won't load `site_libs`), try each widget and anchor link, and confirm
  every figure in `Week<N>_<Topic>_files/` is tracked (`git add -f`).

## Companion R script

`lecture_notes/CODE_Week<N>_<Topic>.R` (see
`CODE_Week10_ModelComparison.R`), linked from the page and posted to
WebCampus.

- The page's visible code, in page order, with one subheading per page
  section (`# 4. AIC ----`) under a header per H1. One- or two-line
  title block (title, packages).
- Brief comments per major step; equations where they say it better.
  `#  [PAUSE] <question>` after the code a page question depends on,
  never the answer.
- Same seeds as the page; run it end to end and check its numbers match
  the page.

## Instructor notes, one per class meeting

`lecture_notes/NOTES_Week<N>_Lecture<M>_<Topic>.Rmd` (Monday = 1,
Wednesday = 2) and `NOTES_Lab<N>_<Topic>.Rmd` (the lab's lecture
portion plus an Act-by-Act quick reference). PDF, `fontsize: 12pt`,
raw LaTeX fine. Model: `NOTES_Week10_WebPageGuide.Rmd` (written as one
doc for the whole week; from Week 11 on, split by meeting).

- **Short and bulleted:** delivery tips that follow the page's
  sections, with time budgets summing to ~30-35 min; what to skip if
  behind.
- **Every page question and discussion prompt**, quoted in full, with a
  one-line *Land:*. Label numbers with what they are ("lasso says 0.07;
  truth 0.3").
- **Widget walkthroughs:** the slider stops to use, in order, the
  values each stop shows, and where to ask students to predict before
  moving it.
- Board work the page doesn't show, and what not to give away before
  lab.
- No "How to prepare" block (prep goes in Canvas announcements, written
  by the instructor). These PDFs sit in `lecture_notes/` and are
  reachable by URL once pushed; nothing links to them.
- **Reference handouts** are separate: `NOTES_Reference_<Topic>.Rmd`,
  student-facing, writing-style rules apply (see
  `NOTES_Reference_AICPenalty.Rmd`).

## After a lecture ("as taught")

When the instructor reports what was covered: (1) add a dated "As
taught" note to that meeting's instructor notes; (2) move uncovered
items into the next meeting's notes (the page usually stays as is,
unless material moves to another week's page); (3) update
`CourseSchedule.csv` and re-render with
`rmarkdown::render_site("schedule.Rmd")`; (4) update the calibration
paragraph above.

## Backlog: web pages for Weeks 1-9

In the next few weeks, build a web page (and companion script, if the
week lacks one) for each earlier week, with widgets where they help.
Sources: that week's PDF NOTES (Weeks 1-10, kept as the record of what
was planned and taught), CODE scripts, WODs, and legacy `LECTUREn.Rmd`
pages. Ask the instructor before removing legacy `LECTUREn` entries
from the navbar or deleting old NOTES.

# Modeling workflow framing

Use this workflow (adopted Week 8; earlier lectures used a shorter
loop and don't need retrofitting) whenever material presents the
scientific process:

$$\text{PROBLEM} \to \text{THEORIES} \to \text{MODELS} \to \text{FIT}
  \to \text{COMPARE} \to \text{REVISE} \circlearrowleft$$

- PROBLEM (often a management question, written first); THEORIES
  (multiple working hypotheses); MODELS (process + noise for each);
  FIT (likelihood, optimization, Bayes, MCMC); COMPARE; REVISE (back
  to the problem and theories, not just parameters). Dashed side box:
  MODELS -> DATA -> FIT (study design, e.g. power analysis).
- **COMPARE** = (1) goodness of fit, (2) prediction and comparison
  against rivals including the null, (3) relevance to the problem (a
  judgment, not a statistic). Survivors are **corroborated**, not
  proven; failures are **falsified**. AIC weights and averaging rank
  survivors; they don't make one true. Avoid "prove/confirm."
- Each big-picture opening names its box: bestiary, distributions,
  simulation = MODELS; likelihood, optimization, Bayes, MCMC = FIT;
  LRT, AIC, cross-validation, LOO, predictive checks = COMPARE.
- Draw it with `draw_workflow()` (copy from `Week10_ModelComparison.Rmd`;
  `highlight = "COMPARE"`), in a chunk with `dev = "png", dpi = 250`
  (the pdf device mangles hyphens).
- Bolker's Figure 1.5 puts data collection second, but his p. 26
  footnote calls that a concession to having data in hand, so
  problem-first matches his own ideal.

# Notation

- Likelihood: $\mathcal{L}(\theta)$, never
  $\mathcal{L}(\theta \mid \text{data})$ (the bar reads as a
  probability of the parameter, and clashes with the posterior
  $P(\theta \mid y)$). Mark data size with a subscript,
  $\mathcal{L}_n(\theta)$. Genuine conditional probabilities like
  $P(y \mid \theta)$ keep the bar.

# Labs: two Acts

- **Act 1:** closed-laptop math, concepts, or planning, so every group
  finishes one Act in the room. **Act 2:** coding, closed -> closed ->
  open -> closed, and usable as an individual take-home.
- **Time:** ~25 min of opening remarks leaves **~2h20m** for both Acts.
  A paper cycle is ~15 min individual + ~5 min group discussion and
  revision (~20 min); coding Acts run 40+ min.
- **Estimates run short:** multiply gut estimates by ~1.5 for math and
  2+ for coding. Fit the period by cutting parts or moving them to an
  "if you finish early" block, never by shrinking times or space.
- **Calibration:** Lab 5 was about right but only one group finished;
  plan slightly less than Lab 5. Lab 7: most groups got 45+ min on Act
  2 and made good progress; it went home to finish. Update this when
  the instructor reports.
- Weeks with real model-vs-model competition can use a "Model Court"
  debate format instead (not yet used).

**Format** (templates: `labs/Lab5_Act1.Rmd`, `labs/Lab5_Act2.Rmd`):

- YAML: `title: "Lab N, Act X: <Name>"`, a one-line `subtitle` posing
  the driving question, `author: "NRES 746, name: _______________   group: _____________________"`
  (group blank wider: it holds groupmates' names), `date: ""`. Labs 1-4
  use older header forms; don't copy them.
- **Scenario** first: 1-2 paragraphs of ecological framing ending on
  the question; data presented as real measurements with a `kable()`.
- Parts `## Part 1: ...`, items `**1a.**`, each with answer space sized
  to the answer and matched to its kind (lines vs. space).
- **Group discussion and revision** (closing): a one-line prompt to
  compare and revise (optionally naming parts worth comparing), then a
  generic instruction to record changes, insights, and muddiest
  points, and one open `blank_lines(6)` block. No pre-written
  reflection questions.
- **Take-home-capable Act** (usually the last): no group language
  ("Checkpoint" for a check-in, "Reflection" for the discussion), and
  no laptop banners; instead a non-enforced "try this on paper first."
  See `labs/Lab2_Act3.Rmd`.
- Both Acts follow one dataset/narrative, and the last Act's closing
  ties back to Act 1's scenario.

# Answer keys (labs and WODs)

- One per handout, built once the handout is stable:
  `labs/answer_keys/Lab<N>_Act<M>_Key.Rmd`,
  `wod/answer_keys/WOD_Week<N>_Lecture<M>_<Topic>_Key.Rmd` (full topic
  slug). Every WOD gets a key. These differ from `NOTES_Lab` files
  (pre-lab instructor prep).
- **Posted to students** (WOD keys for self-checking; lab keys after
  grading), so the writing-style rules apply. Acts never handed out
  (Lab 4 Act 3) get no posted key.
- **Format:** the handout's full text, order, tables, and figures, with
  answers added where the blanks were (drop the name line and answer
  space). **Answers only**: the answer or a couple of alternatives; a
  one-line sample for open-ended items; no explanations or "watch
  for." Compute every number with executed R code. Exception: a WOD
  ending on a teaser (e.g. the birthday problem) can show the key steps
  of the math, since it's walked through in class.
- **Banner:** `\ans{\textbf{Answer key.} Answers are in magenta
  sans-serif type; longer answers have a bar in the margin.}`
- **Black-and-white safe** (from Oct. 22; no color printer): answers
  are magenta and sans-serif; standalone answers (paragraphs, display
  math, code chunks with `highlight = FALSE`) get a left-margin bar.
  `header-includes`:

  ```
  - \usepackage{xcolor}
  - \usepackage{sansmath}
  - \usepackage{framed}
  - \definecolor{ans}{RGB}{200,0,130}
  - \newcommand{\ans}[1]{\textcolor{ans}{\sffamily\sansmath{}#1}}
  - \newcommand{\ansbegin}{\begin{leftbar}\begingroup\color{ans}\sffamily\sansmath}
  - \newcommand{\ansend}{\endgroup\end{leftbar}}
  ```

  `\ans{}` for inline blanks; `\ansbegin`/`\ansend` on their own lines
  (blank lines around) for everything else. Write `\sansmath{}#1`
  (YAML reads ` #` as a comment), and keep the group inside `leftbar`
  (or the color leaks). In figures and `blank_table()` fills use
  `rgb(200, 0, 130, maxColorValue = 255)` plus a non-color cue (line
  type, symbol, legend label, italic `font = 3` cells). Example:
  `wod/answer_keys/WOD_Week10_Lecture1_ConnectTheDots_Key.Rmd`. Older
  keys are color-only; don't retrofit.

# Difficulty

- Graduate level, but no more math than Calc II and intro stats. Build
  the minimum formal machinery; let ecological intuition carry the
  rest.
- **Minimal scaffolding.** Don't pre-break problems into steps or hint
  at the path; wrong first attempts are fine. Some individual questions
  should stay unresolved until the group phase. When in doubt, cut
  scaffolding.
- Every in-class Act builds in comparing and defending answers, and a
  generic, recorded chance to name muddiest points.

# Before considering a task done

- **Knit** every new `.Rmd` (both formats where feasible) and look at
  the rendered output: page breaks, clipped figures, anticlimactic
  numbers (e.g. a p-value of 0 from too few reps). Exception: during an
  instructor revision pass, edit the `.Rmd` only and knit when asked.
- **Look at figures** before settling on parameter values.
- **Verify numerically** any statistical or mathematical claim that goes
  into student materials.
- **Flag ugly numbers** where students compute by hand; adjust the
  example instead.
- **Plotmath `<=`/`>=`** renders as a stray mark after knitting; put
  inequalities in the prose.
- **"Can't write on file ... .pdf"** is a file lock (Acrobat, or a
  brief indexer lock). Retry; if it persists, render a copy in a scratch
  folder to check content and ask the instructor to close the viewer.
  A failed render can delete the PDF; check it's still there. Delete
  leftover `.log`/`.tex`.
- **LaTeX backslashes:** edit with the Edit tool, not `sed` (it has
  turned `\lambda` into `ambda`).
- **Website HTML:** `.gitignore` excludes `*_files/`, so check every
  figure a new page references is tracked (`git ls-files
  --error-unmatch`); fix with `self_contained: true` or `git add -f`,
  and say which.
- **Renames:** grep the repo for the old filename (`index.Rmd`, this
  file, cross-references between documents).

# Development status (update at the end of each thread)

Snapshot: Tue., Oct. 27, 2026 (Week 10).

- **New format from Week 10:** one web page per week plus short
  per-meeting instructor notes (see "Weekly materials"). Week 10's page
  (`Week10_ModelComparison.Rmd`, in the navbar) and delivery guide
  (`NOTES_Week10_WebPageGuide.Rmd`) are built; the Week 10 PDF NOTES
  are now the planning record. Backlog: pages for Weeks 1-9.

- **Weeks 1-9 built and taught.** Week 7 and Week 8 Wednesday NOTES
  carry "as taught" notes; Week 8 Monday has no report. Every WOD and
  every handed-out lab Act through Lab 7 has a (color-only) key.
- **Lab 7 data** (Lab 8 builds on it): Act 1 is one pond, 2 of 20
  frogs Bd-positive, Beta(2, 8) prior; Act 2 is a 12-pond survey
  (`data/BdPondSurvey.csv`, simulated but presented as field data; P11
  is Act 1's pond), logistic infection vs. water temperature, with
  centering at 17 °C fixing a badly mixing sampler.
- **Week 9 as taught:** Lab 7 Act 2 went home after 45+ min of group
  work; the HMC/Stan lecture went well (it used the myxo Gamma model,
  not chytrid). Some students still lack a
  working Stan install (take-home WOD due Mon., Oct. 26); fix installs
  at the start of Lab 8.
- **Week 10 (Oct. 26-28): model comparison and regularization.** Theme:
  every method is a penalty (AIC on the number of parameters; ridge and
  lasso on coefficient size; shrinkage priors as the Bayesian version;
  CV and LOO on held-out data). Reading: Bolker 6.6 (printed p.
  267-290). Tredennick et al. 2021 (Ecology 102: e03336) and Hooten &
  Hobbs 2015 (Ecol. Monogr. 85: 3-28) were optional, not assigned.
  - **Built (Monday):** NOTES, `CODE_Week10_ModelComparison.R`
    (Monday), WOD *Connect the Dots* + key. Order: big picture, WOD,
    lecture. Fir models for LRT/AIC (no-wave best; full vs. no-wave
    D = 2.01, p = 0.57); simulated brown creepers (50 stations, 10
    covariates, truth = big trees, snags, road distance; seed 19) for
    all-subsets AIC vs. ridge/lasso (`glmnet`). The notes correct
    legacy LECTURE8: Poisson *is* nested in the negative binomial (on
    the boundary).
  - **Built (Lab 8):** `labs/Lab8_Act1.Rmd`, `Lab8_Act2.Rmd`, both
    keys, `NOTES_Lab8_ModelComparisonAndRegularization.Rmd`, data
    `data/BdRegionalSurvey.csv` (60 ponds, 8 standardized features;
    simulated, presented as field data; truth = temp, bullfrogs,
    canopy; temp-elevation r = -0.94; generator `sim_ponds(3)` is not
    in the repo). Act 1 (paper): AIC table for 7 models (Model 5 =
    truth is best; Model 6 within 2 but an uninformative parameter),
    two LRTs, temp vs. elevation inseparable, and ridge/lasso budgets
    at t = 0.3 on temp-canopy contours (ridge about (-0.26, 0.14);
    the lasso hits the corner and zeros canopy; 2b asks whether ridge
    can ever zero a coefficient, 2c when the lasso does; changed from a
    temp-elevation panel on Oct. 22). Act 2 (take-home capable): ridge NLL by hand with
    `optim()` (temp -0.64 at lambda 0, -0.41 at 10), ridge path,
    `cv.glmnet` lasso (keeps temp, elev, canopy, bullfrogs), CV
    pseudocode; no Stan. Tuesday opens with a Stan check (~10 min) and
    a ~34-min pre-lab lecture; its catch-up section is a placeholder
    to fill in once Monday's coverage is known. Push the CSV to
    `gh-pages` before the handout goes out.
  - **Built (Wednesday), before Monday was taught; adjust after
    Monday's report:** `NOTES_Week10_Lecture2_BayesianModelComparison.Rmd`,
    the script's Wednesday section (Stan; run Monday's section first),
    no in-class WOD (the ridge-as-prior derivation, lambda = 1/(2
    tau^2), lasso = Laplace prior, is a ~6-min board + R demo in
    section 2; it was a WOD until Oct. 22), and take-home WOD *Comparing Pond Models with LOO*
    + key (students fill in the `log_lik` line; four pond models; B =
    temp + bullfrogs + canopy is best, looic 221.8; D = all eight with
    a Normal(0, 0.3) prior is 2.8 elpd behind, SE 2.7, with p_loo 8.2
    vs. 9.7 for the vague prior). Lecture: fir no-wave vs. full by AIC,
    DIC, WAIC, LOOIC (2278.0/2278.0/2278.6/2278.6 vs. about 2282; elpd
    difference 2.1, SE 1.8); creeper slopes under vague, Normal(0,
    0.2), and Laplace(0, 0.15) priors (northness 0.23 -> 0.08);
    horseshoe as a sketch only; Bayes factors and RJMCMC as a bonus section.
    Posterior predictive checks move to Week 11.
- **Unused material, free to reuse:** Lab 4 Act 3 (raw vs. log-scale
  analysts); Bolker's 7.4 fitting checklist. Lab 4 Act 2 was only a
  live-coding demo, without its AI part (the instructor wants an AI
  live-coding demo in a later lecture). *Curvature by Hand*
  (`WOD_Week8_Lecture2_CurvatureByHand.Rmd`) went out as optional and
  uncollected, so don't assume students did it.
- **Framing to keep:** Hessian intervals are the practical default for
  multiparameter models; profiles are the check (Bolker's cautions
  apply). Gamma regressions use Bolker's parameterization (constant
  shape, scale = mean/shape).
