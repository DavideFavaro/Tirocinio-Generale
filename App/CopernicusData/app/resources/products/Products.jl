module Products

import SearchLight: AbstractModel, DbId
import Base: @kwdef

export Product

@kwdef mutable struct Product <: AbstractModel
  id::DbId = DbId()
end

end
