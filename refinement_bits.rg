--Copyright (c) 2018, Triad National Security, LLC
--All rights reserved.

--This program was produced under U.S. Government contract 89233218CNA000001 for
--Los Alamos National Laboratory (LANL), which is operated by Triad National
--Security, LLC for the U.S. Department of Energy/National Nuclear Security
--Administration.

--THIS SOFTWARE IS PROVIDED BY TRIAD NATIONAL SECURITY, LLC AND CONTRIBUTORS "AS
--IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
--IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
--DISCLAIMED. IN NO EVENT SHALL TRIAD NATIONAL SECURITY, LLC OR CONTRIBUTORS BE
--LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
--CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
--SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
--INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
--CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
--ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
--POSSIBILITY OF SUCH DAMAGE.

--If software is modified to produce derivative works, such modified software should be
--clearly marked, so as not to confuse it with the version available from LANL.
import "regent"

-- keeps track of refinement state of the binary tree
-- might be more performant if a bitfield was used
fspace RefinementBits
{
  isActive: bool,
  isRefined: bool,
  needsRefinement: bool,
  cascadeRefinement: bool,
  wantsCoarsening: bool,
  plusXMoreRefined: bool,
  minusXMoreRefined: bool,
  plusXMoreCoarse: bool,
  minusXMoreCoarse: bool
}

