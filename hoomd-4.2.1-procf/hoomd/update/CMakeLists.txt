### Modified by PRO-CF #[PROCF2023] ###

set(files __init__.py
          box_resize.py
          box_shear.py #[PROCF2023]
          remove_drift.py
          custom_updater.py
          particle_filter.py
   )

install(FILES ${files}
        DESTINATION ${PYTHON_SITE_INSTALL_DIR}/update
       )

copy_files_to_build("${files}" "update" "*.py")
