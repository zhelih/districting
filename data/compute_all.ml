(* from ExtLib *)
let split str sep =
  let p = String.index str sep in
  let len = 1 in
  let slen = String.length str in
  String.sub str 0 p, String.sub str (p + len) (slen - p - len)

let nsplit str sep =
  if str = "" then []
  else
    let rec loop acc pos =
      if pos > String.length str then
        List.rev acc
      else
        let i = try String.index_from str pos sep with Not_found -> String.length str in
        loop (String.sub str pos (i - pos) :: acc) (i + 1)
    in
    loop [] 0

open Printf
let () =
  (* compute graph from <state>_features.csv *)
  if Array.length Sys.argv < 2 then begin
    printf "Usage: %s <state>\n" Sys.argv.(0);
    exit 0
  end;
  let state = Sys.argv.(1) in

  let f_in = try open_in (state ^ "_features.csv") with _ -> (printf "Input not found\n"; exit 1) in
  (* format: geoid, pop, nbs *)
  (* skip first row *)
(*   ignore (input_line f_in); *)

  (* can be done in one read, but I am toooo lazy *)
  let h_pop = Hashtbl.create 1 in
  let h_ind = Hashtbl.create 1 in
  let rec loop i =
    try
      let line = input_line f_in in
(*       printf "1: Processing line %S\n" line; *)
      let (geoid,rest) = split line ',' in
      let (pop,rest) = split rest ',' in
      let pop = try int_of_string pop with _ -> 0 in
      Hashtbl.add h_pop geoid pop;
      Hashtbl.add h_ind geoid i;
      loop (i+1)
    with End_of_file -> ()
  in
  loop 1;
  seek_in f_in 0;

  (* load distance matrix (memory intence) *)
  let f = open_in (state ^ "_d.csv") in
  let header = input_line f in
  let fields = match nsplit header ',' with "ID"::tl -> tl | _ -> raise (Failure "Bad distance CSV") in
  let nr_col = List.length fields in
  let h_col = Hashtbl.create nr_col in
  let h_row = Hashtbl.create nr_col in
  List.iteri (fun i s ->
    Hashtbl.add h_col (String.trim s) i
  ) fields; (* maps geoid to matrix index *)

  let distances = Array.make_matrix (nr_col+1) (nr_col+1) 0. in
  let rec loop i =
    try
      let (geoid,d) = match nsplit (input_line f) ',' with f::s -> f,s | _ -> assert false in

      Hashtbl.add h_row geoid i;
      List.iteri (fun j d ->
        distances.(j).(i) <- float_of_string (String.trim d)
      ) d;
      loop (i+1)
    with End_of_file -> ();
  in
  loop 0;
  close_in f;

(*   printf "first loop done\n"; *)

  let n = Hashtbl.length h_ind in (* nr_nodes *)
  let adj = Array.make_matrix (n+1) (n+1) None in
  (* skip first row *)
  ignore (input_line f_in);
  let rec loop i =
    try
      let line = input_line f_in in
(*       printf "2: Processing line %S\n" line; *)
      let (geoid,rest) = split line ',' in (* too lazy to merge with code above *)
      let (_,rest) = split rest ',' in
      let rest = if rest.[0] = '"' then String.sub rest 1 (String.length rest - 2) else rest in (* remove quotes *)
      let nbs = nsplit rest ',' in

      let v_from = try Hashtbl.find h_ind geoid with exn -> (printf "missing %s (v_from)\n" geoid; raise exn) in
      List.iter (fun geoid_to ->
        let v_to = try Hashtbl.find h_ind geoid_to with exn -> (printf "missing %s from %s (v_to)\n" geoid_to geoid; raise exn) in
        let d_col = try Hashtbl.find h_col geoid with exn -> (printf "missing %s from distances cols" geoid; raise exn) in
        let d_row = try Hashtbl.find h_row geoid_to with exn -> (printf "missing %s from distances rows" geoid_to; raise exn) in

        (* some checks *)
        let d_col2 = try Hashtbl.find h_col geoid_to with exn -> (printf "missing %s from distances cols" geoid_to; raise exn) in
        let d_row2 = try Hashtbl.find h_row geoid with exn -> (printf "missing %s from distances rows" geoid; raise exn) in
        let d = distances.(d_col).(d_row) in
        let d2 = distances.(d_col2).(d_row2) in
        if d <> d2 then
          printf "Distance error: %.8f <> %.8f for geoid %s and %s (indices (%d,%d) and (%d,%d)\n" d d2 geoid geoid_to d_col d_row d_col2 d_row2;
        adj.(v_from).(v_to) <- Some d;
        adj.(v_to).(v_from) <- Some d
      ) nbs;
      loop (i+1)
    with End_of_file -> ()
  in
  loop 1;
  close_in f_in;

  let f = open_out (state ^ ".dimacs") in

  fprintf f "p edge %d\n" n;
  for i = 1 to n do
    for j = i+1 to n do
      match adj.(i).(j) with
(*       | Some d -> fprintf f "e %d %d %.8f\n" (i-1) (j-1) d *)
      | Some d -> fprintf f "e %d %d\n" (i-1) (j-1)
      | None -> ()
    done
  done;
  close_out f;

  printf "adj done\n";

  let f = open_out (state ^ ".hash") in
  Hashtbl.iter (fun k v ->
    fprintf f "%d %s\n" (v-1) k
  ) h_ind;
  close_out f;
  let f = open_out (state ^ ".population") in
  (* compute total population *)
  let total_pop = ref 0 in
  Hashtbl.iter (fun _k v -> total_pop := !total_pop + v) h_pop;
  fprintf f "total pop = %d\n" !total_pop;
  Hashtbl.iter (fun k v ->
    let ind = Hashtbl.find h_ind k in
    fprintf f "%d %d\n" (ind-1) v
  ) h_pop;
  close_out f
