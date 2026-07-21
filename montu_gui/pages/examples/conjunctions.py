# %pip install montu
"""Planetary conjunction search with MontuPython."""

import montu

mars = montu.Planet("Mars")
aldebaran = montu.Stars(subset="bright", ProperName="Aldebaran", return_as="Star")

explorer = montu.ConjunctionExplorer(
    bodies=[mars, aldebaran],
    maxseparation=5,
)

start = montu.Time("2022-09-01")
end = montu.Time("2022-10-01")
hits = explorer.search(start=start, end=end, observer="geocentric", verbose=False)

print(f"Conjunctions found: {len(hits)}")
for conj in hits:
    print(
        f"  {conj.mtime.readable.datemix}  "
        f"sep={conj.separation:.2f}°  "
        f"bodies={'–'.join(conj.body_names)}"
    )
