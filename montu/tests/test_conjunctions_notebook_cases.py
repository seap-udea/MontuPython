"""Ground-truth conjunction cases from examples/MontuPython-Conjunctions.ipynb.

Tolerances (as requested):
- angular separation at closest approach: 0.1°
- epoch of closest approach: 3 days
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Sequence

import pytest

import montu

SEP_TOL_DEG = 0.1
DATE_TOL_DAYS = 3.0


def _star(proper_name: str):
    return montu.Stars(subset="bright", ProperName=proper_name, return_as="Star")


def _parse_time(text: str, calendar: str | None = None) -> montu.Time:
    if calendar:
        return montu.Time(text, calendar=calendar)
    return montu.Time(text)


def _search(
    bodies: Sequence,
    maxseparation: float,
    start: str,
    end: str,
    *,
    calendar: str | None = None,
) -> list[montu.Conjunction]:
    start_t = _parse_time(start, calendar)
    end_t = _parse_time(end, calendar)
    explorer = montu.ConjunctionExplorer(bodies=bodies, maxseparation=maxseparation)
    return explorer.search(
        start=start_t,
        end=end_t,
        observer="geocentric",
        verbose=False,
    )


def _assert_hit_near(
    hit: montu.Conjunction,
    expected_date: str,
    expected_separation: float,
    *,
    calendar: str | None = None,
) -> None:
    expected = _parse_time(expected_date, calendar)
    delta_days = abs(hit.mtime.jed - expected.jed)
    assert delta_days <= DATE_TOL_DAYS, (
        f"epoch delta {delta_days:.2f} d > {DATE_TOL_DAYS} d "
        f"(got {hit.mtime.readable.datespice}, expected {expected.readable.datespice})"
    )
    assert hit.separation == pytest.approx(expected_separation, abs=SEP_TOL_DEG)
    assert hit.in_conjunction


def _match_expected_hits(
    hits: Sequence[montu.Conjunction],
    expected: Sequence[tuple[str, float]],
    *,
    calendar: str | None = None,
) -> None:
    assert len(hits) == len(expected)
    used: set[int] = set()
    for expected_date, expected_sep in expected:
        expected_jed = _parse_time(expected_date, calendar).jed
        index, best = min(
            ((idx, hit) for idx, hit in enumerate(hits) if idx not in used),
            key=lambda pair: abs(pair[1].mtime.jed - expected_jed),
        )
        used.add(index)
        _assert_hit_near(best, expected_date, expected_sep, calendar=calendar)


@dataclass(frozen=True)
class NotebookConjunctionCase:
    category: str
    bodies_factory: Callable[[], Sequence]
    maxseparation: float
    start: str
    end: str
    calendar: str | None
    expected: Sequence[tuple[str, float]]


NOTEBOOK_CASES: tuple[NotebookConjunctionCase, ...] = (
    NotebookConjunctionCase(
        category="planet_star",
        bodies_factory=lambda: [montu.Planet("Mars"), _star("Aldebaran")],
        maxseparation=5.0,
        start="2022-09-01",
        end="2022-10-01",
        calendar=None,
        expected=(("2022-09-07", 4.2746),),
    ),
    NotebookConjunctionCase(
        category="planet_planet",
        bodies_factory=lambda: [montu.Planet("Jupiter"), montu.Planet("Saturn")],
        maxseparation=5.0,
        start="bce 7-01-01",
        end="bce 7-12-31",
        calendar="mixed",
        expected=(
            ("bce 7-05-27", 0.9846),
            ("bce 7-09-28", 0.9745),
            ("bce 7-12-03", 1.0541),
        ),
    ),
    NotebookConjunctionCase(
        category="planet_planet_planet",
        bodies_factory=lambda: [
            montu.Planet("Mercury"),
            montu.Planet("Mars"),
            montu.Planet("Saturn"),
        ],
        maxseparation=5.0,
        start="2026-04-15",
        end="2026-04-25",
        calendar=None,
        expected=(("2026-04-20", 1.6511),),
    ),
    NotebookConjunctionCase(
        category="planet_planet_planet",
        bodies_factory=lambda: [
            montu.Planet("Venus"),
            montu.Planet("Jupiter"),
            montu.Planet("Mercury"),
        ],
        maxseparation=5.0,
        start="2013-05-20",
        end="2013-05-31",
        calendar=None,
        expected=(("2013-05-27", 2.3634),),
    ),
    NotebookConjunctionCase(
        category="planet_planet_star",
        bodies_factory=lambda: [
            montu.Planet("Venus"),
            montu.Planet("Mars"),
            _star("Regulus"),
        ],
        maxseparation=5.5,
        start="2021-07-15",
        end="2021-07-25",
        calendar=None,
        expected=(("2021-07-22", 4.9933),),
    ),
)


@pytest.mark.parametrize(
    "case",
    NOTEBOOK_CASES,
    ids=[
        "mars-aldebaran-2022",
        "jupiter-saturn-7-bce",
        "mercury-mars-saturn-2026",
        "venus-jupiter-mercury-2013",
        "venus-mars-regulus-2021",
    ],
)
def test_notebook_conjunction_cases(case: NotebookConjunctionCase):
    hits = _search(
        case.bodies_factory(),
        case.maxseparation,
        case.start,
        case.end,
        calendar=case.calendar,
    )
    _match_expected_hits(hits, case.expected, calendar=case.calendar)
