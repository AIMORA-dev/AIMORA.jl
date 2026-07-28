module ArcFlashStudy

using ..StudyCore

export ARC_FLASH_STUDY, run_arc_flash

const ARC_FLASH_STUDY = StudyDescriptor(
    :arc_flash,
    "Arc flash",
    :safety,
    :planned,
    "src/studies/arc_flash.jl",
)

run_arc_flash(args...; kwargs...) = study_not_implemented(ARC_FLASH_STUDY)

end
