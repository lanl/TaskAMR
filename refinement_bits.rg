import "regent"

-- keeps track of refinement state of the binary tree
-- might be more performant if a bitfield was used
fspace RefinementBits
{
  isActive: bool,
  isRefined: bool,
  needsRefinement: bool,
  wantsCoarsening: bool,
  plusXMoreRefined: bool,
  minusXMoreRefined: bool,
  plusXMoreCoarse: bool,
  minusXMoreCoarse: bool
}

