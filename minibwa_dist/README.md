# minibwa_dist — maintaining the distribution

`nh13/minibwa` is a **downstream distribution** of [lh3/minibwa](https://github.com/lh3/minibwa):
upstream, plus the changes upstream declined. This directory holds the machinery that keeps it
that way. `README.md` at the repo root is for people who *use* the distribution; this file is for
whoever *maintains* it.

## The model, in one paragraph

`minibwa_dist/features.toml` is the single source of truth. It lists feature branches. A sync run
takes current upstream `master`, **merges** each listed branch onto it in order, and the result is
`dist`. Nothing is ever merged *into* `dist` — it is rebuilt from scratch and force-pushed. Adding
a feature means adding a manifest block; removing one means deleting that block.

```
lh3/master ──► master (mirror) ──┐
                                 ├─ merge each [[feature]] in order ─► dist-next ──ship──► dist ──► tag
nh13/feat/*  (frozen branches) ──┘
```

**Merge, never rebase.** A feature branch stays frozen at the commit it branched from, so its side
of any conflict is byte-stable and `git rerere` can replay the resolution on every later run.
Rebasing a feature to "keep it current" shifts that content and defeats the cache — it is the one
change that quietly makes this system expensive.

## Branch rules

| branch | what it is | may I commit to it? |
|---|---|---|
| `master` | pure mirror of `lh3/master`, fast-forwarded by the bot | **no** — merging into it breaks the assembly base |
| `feat/*`, `perf/*`, `chore/*`, `fix/*` | frozen feature branches | yes, that is where features live |
| `ops/distro` | this tooling, the manifest, the workflows | yes |
| `dist` | rebuilt artifact, force-pushed every sync | **no** — anything here is destroyed on the next run |
| `dist-next` | the candidate under review | **no** |
| `rerere-cache` | orphan branch of conflict resolutions | only via the tooling |

Consume **release tags**, not branches. `dist` rewinds.

## Adding a feature

```sh
# 1. Branch off UPSTREAM master -- never off dist. The merge base must stay fixed.
cd ~/work/git/minibwa/main
git fetch origin master
git worktree add ../my-feature -b feat/my-thing origin/master
cd ../my-feature
# ... write it, commit ...

# 2. Push to the fork. The manifest resolves nh13/<branch>.
git push nh13 feat/my-thing
```

3. Add a block to `minibwa_dist/features.toml` on `ops/distro`:

```toml
[[feature]]
name     = "my-thing"
branch   = "feat/my-thing"
required = false                 # only the tooling entry is ever true
output   = "identical"           # see the table below
summary  = "what it does, and the measured win if there is one"
upstream = { status = "unsubmitted" }
```

Merge order is the manifest order, and it is deliberate: features that touch no shared file go
first, the hub feature (`alt-liftgroup`, which touches nine) goes last. Put a new feature where its
file overlap suggests, not at the end by default.

`required = false` means a merge conflict **drops** the feature from that build and files an issue,
rather than failing the build. Only `ops-distro` is `required = true`, because without it the
workflows are not on `dist` at all.

## The `output` classes, and what CI enforces

| value | meaning | what the gate does |
|---|---|---|
| `identical` | byte-identical to stock upstream, always | aligns the chrM fixture with stock and with the build; **any** differing SAM record fails |
| `conditional` | identical unless a named condition holds | requires `condition` and `negative`; CI runs the **negative** case and demands byte-identity |
| `changes-output` | unconditionally different | **blocked from release** until a real-data concordance measurement exists |

```toml
output    = "conditional"
condition = "-L is set"
negative  = "align the chrM fixture without -L; SAM must be byte-identical to stock"
```

The gate is the only thing standing between a replayed-but-now-wrong conflict resolution and a
silently different aligner. Declaring `identical` and being wrong will fail the build, which is the
intended outcome.

Note what is **not** covered: the *positive* path of a `conditional` feature. Those change output
on purpose, so byte-identity says nothing about them; they rely on each feature branch's own tests.

## `upstream.status` — six values, each with consequences

| status | meaning | what the automation does |
|---|---|---|
| `unsubmitted` | never offered upstream | listed in a monthly "candidates to offer" issue |
| `open` | PR open upstream | carried; the branch is never auto-rebased (that would invalidate the PR) |
| `rejected` | upstream declined | permanent resident; never nominated |
| `merged` / `superseded` | upstream now carries it | **rejected by the parser** — the block must be deleted |
| `withdrawn` | investigated and disproven | lives in `[[withdrawn]]`, never in the build, never nominated |

`withdrawn` exists so a disproven idea stays disproven. Each entry carries a `report` path; see
`GRAVEYARD.md`. Do not "helpfully" resurrect one without reading its report first.

## Removing a feature

Delete its `[[feature]]` block. That is the whole operation.

If upstream *takes* a feature, the weekly reconciler notices the PR went merged and opens a PR
deleting the block for you. Merge it. Leaving a graduated feature in the manifest guarantees a
duplicate-application conflict on every run.

## When a sync drops a feature

The bot files an issue naming the feature and the conflicting files, ships the rest, and carries
on. To fix it, train `rerere` **at the exact intermediate state a later run reproduces** — which is
what `--stop-at` gives you:

```sh
cd ~/work/git/minibwa/ops-distro
git fetch nh13 '+refs/heads/*:refs/remotes/nh13/*' --prune
git config rerere.enabled true && git config rerere.autoUpdate true

pixi exec -s python=3.12 python -m minibwa_dist.cli sync --repo . --base master \
  --out dist-next --remote nh13 --stop-at <feature>     # leaves the conflict in the worktree

git diff --name-only --diff-filter=U     # resolve every one of these
make -j8                                  # MUST build before you record anything
git add -- <resolved paths>               # explicit paths only
git -c commit.gpgsign=false commit --no-edit
git reset --hard HEAD~1                   # discard the merge; rerere keeps the resolution
```

You should see `Recorded resolution for '<file>'`. Re-run the plain sync; the feature should now
merge. Then publish the cache:

```sh
git push nh13 rerere-cache
```

Resolve **in manifest order** when several features are dropped: each one you fix changes the
"ours" side for the ones after it.

**How to resolve.** These are usually not semantic disagreements. A feature branch was written
against an older upstream, and the file has moved since — so the answer is normally "take
upstream's evolved code and re-apply the feature's change on top of it". Keep both sides' intent. A
feature that merges but no longer does its job is worse than a dropped one.

## The rerere cache

It lives on the orphan branch `rerere-cache` because `.git/rr-cache` is local and never pushed — a
fresh runner would start cold and re-hit every conflict.

**It never evicts.** `save()` only ever adds. So a resolution that is textually valid but
semantically wrong will replay on every subsequent run and fail the gate every time, with no code
path to detect or expire it. The remedy is manual:

```sh
git fetch nh13 rerere-cache:rerere-cache
git checkout rerere-cache
git rm -r <offending-hash-dir> && git commit -S -m "chore(rerere): drop a bad resolution"
git push --force-with-lease nh13 rerere-cache
```

Find the offending directory by resolving the conflict correctly by hand and comparing.

## Releases

Upstream publishes `v0.6`; we publish `v0.6-nh13.1`. Re-cutting against the same upstream version
increments the revision: `v0.6-nh13.2`.

```sh
gh workflow run distro-release.yml --repo nh13/minibwa -f upstream_tag=v0.6 -f revision=1
```

The job builds, gates and verifies the version string **before** it creates the tag — a tag on an
unverified commit is worse than no release. Release notes are upstream's, plus the feature table,
plus any feature dropped from that build (read back from `dist-manifest.json`, which every `dist`
commit carries).

`MB_VERSION` carries the downstream identity into the SAM `@PG VN:` tag, so any output file is
traceable to the build that produced it. Sync stamps `<upstream>-nh13.dev+<sha>`; release
re-stamps the final version.

## The workflows

| file | trigger | does |
|---|---|---|
| `distro-sync.yml` | daily 06:17 UTC + dispatch | mirror master, assemble `dist-next`, build, gate, open/update the sync PR |
| `distro-ship.yml` | `ship` label on the sync PR | force-update `dist` to the **reviewed** head sha |
| `distro-release.yml` | dispatch | build, gate, verify, then tag and publish |
| `distro-reconcile.yml` | weekly Mon 07:31 UTC | correct stale statuses, PR the manifest, nominate unsubmitted work |
| `distro-test.yml` | push/PR | ruff + the test suite |

### Workflow changes are always one cycle behind

The bot runs the version of itself that is on the **default branch** — `dist` — not the version
on `ops/distro`. GitHub resolves `schedule` and `workflow_dispatch` against the default branch
only. So a change to any `distro-*.yml` does **not** take effect when you push it to
`ops/distro`; it takes effect after the sync that carries it onto `dist` has been shipped.

The practical consequence: **you fix a workflow using the old workflow.** Push the change, let a
sync run (still the old code), ship that PR, and the *next* run uses the new code. Expect the
first run after a workflow change to behave the old way, and do not treat that as the fix having
failed.

The sharper consequence: **a broken workflow cannot repair itself.** If a change lands on `dist`
that breaks `distro-sync`, no future sync can replace it, because syncing is the thing that is
broken. The recovery is to assemble and publish `dist` by hand — the bootstrap sequence — which
is the same reason the first `dist` had to be built manually.

Two habits follow. Keep `distro-test.yml` passing on `ops/distro` before shipping, since it runs
on push there and is the only check that sees a workflow change before `dist` does. And when a
workflow change is risky, ship it on its own rather than bundled with feature changes, so a
manual recovery has a small diff to reason about.

`dist` is advanced in exactly one place — `distro-ship`, on a labelled PR, after the gates passed.
Nothing else writes it.

## Recovering `dist`

`dist` is disposable by construction. If it is wrong, fix the input and re-assemble; there is
nothing to repair in place. Locally: run the sync, verify the gate, then
`git push --force nh13 dist-next:refs/heads/dist`. Existing release tags are unaffected — they
point at their own commits.

## Known limits

Stated plainly so nobody rediscovers them the hard way.

- **rerere is cheaper, not free.** A pre-image is stable only if the *whole prefix* of merge
  outcomes is stable, so the number of distinct pre-images grows with the number of distinct
  drop-subsets, not the number of conflicts. `--stop-at` converges roughly one feature per run.
- **The gate corpus is thin.** Upstream's chrM fixture: ~2,000 records, one contig, all 151 bp
  paired-end, no long reads. `extd2-avx512` is carried on a HiFi/ONT justification and is therefore
  not meaningfully exercised by the gate.
- **The gate runs on x86 only**, while `ll-affine-reassoc` is an arm64-specific kernel change.
- **First run cannot self-start.** Scheduled workflows only fire from the default branch, which is
  `dist`, which does not exist until the first assembly is pushed by hand.
- **`git push origin master:master` in the sync workflow is non-forced.** If upstream ever
  force-pushes `master`, the job wedges until a human intervenes.

## Layout

```
minibwa_dist/
  features.toml     the manifest -- the single source of truth
  manifest.py       parse + validate it (strict on purpose)
  assemble.py       merge features onto a fresh base; stamp the dev version
  rerere_cache.py   persist .git/rr-cache on the orphan branch
  gates.py          byte-identity SAM comparison against stock
  reconcile.py      keep upstream.status honest; delete graduated blocks
  render.py         generate README.md's block, GRAVEYARD.md, release notes
  cli.py            thin entry points; every real decision lives in a module
  gitutil.py        git subprocess wrapper
  tests/            106 tests, no network, synthetic git repos in tmp dirs
```

Python is **stdlib-only** — no runtime dependencies. Run it with
`pixi exec -s python=3.12 -s pytest -s ruff <cmd>`; `tomllib` needs 3.11+.
