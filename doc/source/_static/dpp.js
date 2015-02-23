function run_dpp_demo(size_multiplier,
        font_size_multiplier,
        num_variables,
        image_path) {
    var nv = get_num_variables(num_variables);
    var image_width = 1358;
    var image_height = 1388;
    if (nv == 3) {
        var image_width = 600;
        var image_height = 450;
    }
    size_multiplier = typeof size_multiplier !== 'undefined' ?
            size_multiplier : 0.7;
    font_size_multiplier = typeof font_size_multiplier !== 'undefined' ?
            font_size_multiplier : 1.0;
    create_dpp_form(nv, image_path, font_size_multiplier);
    var a = validate_dpp_form(nv);
    var medium = create_dpp_canvas(
            image_width*size_multiplier, 
            image_height*size_multiplier,
            nv,
            image_path);
    update_dpp_tree(nv, image_path, font_size_multiplier);
}

function create_dpp_form(num_variables, image_path, font_size_multiplier) {
    var nv = get_num_variables(num_variables);
    var dpp_div = document.createElement("div");
    dpp_div.setAttribute("name", "dpp_" + nv + "_div");
    document.getElementsByTagName("body")[0].appendChild(dpp_div);

    dpp_form = document.createElement("form");
    dpp_form.setAttribute("name", "dpp_" + nv + "_form");
    l = document.createElement("label");
    l.for = "cparam";
    l.innerHTML = "Concentration parameter: ";
    a = document.createElement("input");
    a.setAttribute("label", "label");
    a.setAttribute("type", "number");
    a.setAttribute("name", "concentration_param");
    a.setAttribute("id", "cparam");
    a.setAttribute("min", "0.0");
    a.setAttribute("max", "100000.0");
    a.setAttribute("step", "any");
    a.setAttribute("value", "1.5");
    a.setAttribute("onkeypress", "parse_key_press(event, " + nv + ");");
    b = document.createElement("input");
    b.setAttribute("type", "button");
    b.setAttribute("id", "update_" + nv + "_button");
    b.setAttribute("value", "Update");
    b.setAttribute("onclick", "update_dpp_tree(" + nv + 
            ",'" + image_path + "'" +
            "," + font_size_multiplier + ");");
    i = document.createElement("input");
    i.setAttribute("type", "text");
    i.setAttribute("name", "StackOverflow1370021");
    i.setAttribute("value", "Fix IE bug");
    i.setAttribute("style", "display: none;");

    dpp_form.appendChild(l);
    dpp_form.appendChild(a);
    dpp_form.appendChild(b);
    dpp_form.appendChild(i);
    dpp_div.appendChild(dpp_form);
}

function parse_key_press(event, num_variables) {
    var nv = get_num_variables(num_variables);
    if (typeof event == 'undefined' && window.event) {
        event = window.event;
    }
    if (event.keyCode == 13) {
        document.getElementById("update_" + nv + "_button").click();
    }
    return false;
}

function get_num_variables(num_variables) {
    if (num_variables != 3) {
        num_variables = 4;
    }
    return num_variables;
}

function validate_dpp_form(num_variables) {
    var nv = get_num_variables(num_variables);
    var a = parseFloat(document.forms["dpp_" + nv +
            "_form"]["concentration_param"].value);
    if (a == null || a == "" || a < 0.0 || a > 10000000.0) {
        alert("Concentration parameter must be between 0.0 and 10000000.0");
        return false;
    }
    return a;
}

function update_dpp_tree(num_variables, image_path, font_size_multiplier) {
    var nv = get_num_variables(num_variables);
    var a = validate_dpp_form(nv);
    if (nv == 3) {
        annotate_dpp_3_tree(a, image_path, font_size_multiplier);
    } else {
        annotate_dpp_4_tree(a, image_path, font_size_multiplier);
    }
    return false;
}

function create_dpp_canvas(w, h, num_variables, image_path) {
    var nv = get_num_variables(num_variables);
    var canvas = document.createElement("canvas");
    canvas.id = "dpp_" + nv + "_canvas";
    canvas.width = w;
    canvas.height = h;
    canvas.style.border = "1px solid #d3d3d3";
    document.body.appendChild(canvas);
    var context = canvas.getContext("2d");
    embed_dpp_tree_in_canvas(canvas, context, image_path);
    return [canvas, context];
}

function annotate_dpp_3_tree(a, image_path, font_size_multiplier) {
    var canvas = document.getElementById("dpp_3_canvas");
    var context = canvas.getContext("2d");
    var w = canvas.width;
    var h = canvas.height;
    reset_dpp_canvas(canvas, context, image_path);
    context.font = 24*font_size_multiplier + "px Arial";
    context.textAlign="left";
    context.fillText(a, 0.56*w, 0.045*h);
    context.font = 22*font_size_multiplier + "px Arial";
    context.fillText(calc_dpp_prob(a, 1, 2).toFixed(3), 0.905*w, 0.139*h);
    context.fillText(calc_dpp_prob(a, 1, a).toFixed(3), 0.905*w, 0.339*h);
    context.fillText(calc_dpp_prob(a, a, 1).toFixed(3), 0.905*w, 0.559*h);
    context.fillText(calc_dpp_prob(a, a, 1).toFixed(3), 0.905*w, 0.759*h);
    context.fillText(calc_dpp_prob(a, a, a).toFixed(3), 0.905*w, 0.959*h);
    return false;
}

function annotate_dpp_4_tree(a, image_path, font_size_multiplier) {
    var canvas = document.getElementById("dpp_4_canvas");
    var context = canvas.getContext("2d");
    var w = canvas.width;
    var h = canvas.height;
    reset_dpp_canvas(canvas, context, image_path);
    context.font = 24*font_size_multiplier + "px Arial";
    context.textAlign="left";
    context.fillText(a, 0.53*w, 0.02*h);
    context.font = 22*font_size_multiplier + "px Arial";
    context.fillText(calc_dpp_prob(a, 1, 2, 3).toFixed(3), 0.92*w, 0.067*h);
    context.fillText(calc_dpp_prob(a, 1, 2, a).toFixed(3), 0.92*w, 0.133*h);
    context.fillText(calc_dpp_prob(a, 1, a, 2).toFixed(3), 0.92*w, 0.198*h);
    context.fillText(calc_dpp_prob(a, 1, a, 1).toFixed(3), 0.92*w, 0.263*h);
    context.fillText(calc_dpp_prob(a, 1, a, a).toFixed(3), 0.92*w, 0.328*h);
    context.fillText(calc_dpp_prob(a, a, 1, 1).toFixed(3), 0.92*w, 0.393*h);
    context.fillText(calc_dpp_prob(a, a, 1, 2).toFixed(3), 0.92*w, 0.459*h);
    context.fillText(calc_dpp_prob(a, a, 1, a).toFixed(3), 0.92*w, 0.524*h);
    context.fillText(calc_dpp_prob(a, a, 1, 2).toFixed(3), 0.92*w, 0.590*h);
    context.fillText(calc_dpp_prob(a, a, 1, 1).toFixed(3), 0.92*w, 0.655*h);
    context.fillText(calc_dpp_prob(a, a, 1, a).toFixed(3), 0.92*w, 0.722*h);
    context.fillText(calc_dpp_prob(a, a, a, 1).toFixed(3), 0.92*w, 0.787*h);
    context.fillText(calc_dpp_prob(a, a, a, 1).toFixed(3), 0.92*w, 0.853*h);
    context.fillText(calc_dpp_prob(a, a, a, 1).toFixed(3), 0.92*w, 0.918*h);
    context.fillText(calc_dpp_prob(a, a, a, a).toFixed(3), 0.92*w, 0.983*h);
    return false;
}

function reset_dpp_canvas(canvas, context, image_path) {
    clear_canvas(canvas, context);
    embed_dpp_tree_in_canvas(canvas, context, image_path);
}

function embed_dpp_tree_in_canvas(canvas, context, image_path) {
    var img = new Image();
    img.src = image_path;
    img.id = "dpp_image_" + image_path
    context.drawImage(img, 0, 0, canvas.width, canvas.height);
}

function clear_canvas(canvas, context) {
    context.save();
    context.setTransform(1, 0, 0, 1, 0, 0);
    context.clearRect(0, 0, canvas.width, canvas.height);
    context.restore();
    context.beginPath();
}

function calc_dpp_prob(alpha, numerator1, numerator2, numerator3) {
    if(typeof(numerator3)==='undefined') {
        return ((numerator1/(alpha+1))*(numerator2/(alpha+2)));
    } else {
        return ((numerator1/(alpha+1))*(numerator2/(alpha+2))*
                (numerator3/(alpha+3)));
    }
}
