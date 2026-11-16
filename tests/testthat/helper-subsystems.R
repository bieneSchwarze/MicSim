# Shared test helpers for subsystems.
#
# The engine takes subsystems already compiled against a codec and a
# registry, which is what add_subsystem() stores in a model. Engine-level
# tests build both by hand, so they compile by hand too.
compiled <- function(subsystem, codec, registry) {
  compile_subsystem(subsystem, codec, registry)
}
