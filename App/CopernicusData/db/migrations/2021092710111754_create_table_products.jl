module CreateTableProducts

import SearchLight.Migrations: create_table, column, primary_key, add_index, drop_table

function up()
  create_table(:products) do
    [
      primary_key()
      column(:column_name, :column_type)
    ]
  end

  add_index(:products, :column_name)
end

function down()
  drop_table(:products)
end

end
