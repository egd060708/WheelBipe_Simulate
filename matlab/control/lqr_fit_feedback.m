function feedback = lqr_fit_feedback(H_array,K_array,title_name,variant_name)
    feedback = fit_feedback(H_array,K_array,3);
    disp(title_name);
    f1 = sprintf(strcat(variant_name,'.a = %d;\n'),feedback(1,1));
    f2 = sprintf(strcat(variant_name,'.b = %d;\n'),feedback(1,2));
    f3 = sprintf(strcat(variant_name,'.c = %d;\n'),feedback(1,3));
    f4 = sprintf(strcat(variant_name,'.d = %d;\n'),feedback(1,4));
    fprintf(f1);
    fprintf(f2);
    fprintf(f3);
    fprintf(f4);
end
