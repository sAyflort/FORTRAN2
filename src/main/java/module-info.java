module ru.bmstu.fortran.fortran2 {
    requires javafx.controls;
    requires javafx.fxml;


    opens ru.bmstu.fortran.fortran2 to javafx.fxml;
    exports ru.bmstu.fortran.fortran2;
}